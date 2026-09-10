"""Gatheral raw SVI volatility smile: parametrization and calibration.

The raw SVI parametrization models total implied variance as a function of
log-moneyness ``k = log(K / F)``:

    w(k) = a + b * ( rho * (k - m) + sqrt( (k - m)^2 + s^2 ) )

where ``w = sigma_BS^2 * t`` is total variance. Implied vol is then
``sigma(k) = sqrt(w(k) / t)``.

Parameters:
    a   : vertical level of variance (a + b*s*sqrt(1-rho^2) >= 0 keeps w >= 0)
    b   : angle / slope of the wings (b >= 0)
    rho : skew, in (-1, 1)
    m   : horizontal shift of the smile minimum
    s   : smoothness / curvature (sigma in Gatheral's notation, s > 0)
"""

import math
from dataclasses import dataclass
from typing import Sequence, Tuple

from .optimize import nelder_mead


@dataclass(frozen=True)
class SVIParams:
    a: float
    b: float
    rho: float
    m: float
    s: float

    def total_variance(self, k: float) -> float:
        """Total implied variance w(k) at log-moneyness k."""
        return self.a + self.b * (self.rho * (k - self.m)
                                   + math.sqrt((k - self.m) ** 2 + self.s ** 2))

    def implied_vol(self, k: float, t: float) -> float:
        """Black-Scholes implied vol at log-moneyness k for expiry t (years)."""
        w = self.total_variance(k)
        if w < 0:
            w = 0.0
        return math.sqrt(w / t)

    def is_arbitrage_free_wings(self) -> bool:
        """Necessary condition: b*(1+|rho|) keeps wings below Lee's slope bound.

        Lee's moment formula caps the large-|k| slope of total variance at 2.
        Raw SVI's asymptotic slopes are b*(1+rho) (right) and b*(1-rho) (left);
        both must be <= 2 to avoid static (butterfly-independent) wing arbitrage.
        """
        return self.b * (1.0 + abs(self.rho)) <= 2.0 + 1e-9


def calibrate_svi(
    ks: Sequence[float],
    total_variances: Sequence[float],
    weights: Sequence[float] = None,
    initial: SVIParams = None,
    max_iter: int = 4000,
) -> Tuple[SVIParams, float]:
    """Fit raw SVI to observed (log-moneyness, total-variance) points.

    Returns ``(params, rmse)`` where rmse is the root-mean-square total-variance
    error. Uses an unconstrained Nelder-Mead over a smooth reparametrization
    that enforces ``b >= 0``, ``s > 0``, and ``rho in (-1, 1)``.
    """
    ks = [float(k) for k in ks]
    tv = [float(v) for v in total_variances]
    n = len(ks)
    if n < 5:
        raise ValueError("SVI needs at least 5 quotes to fit 5 parameters")
    if weights is None:
        weights = [1.0] * n
    wsum = sum(weights)

    # Reasonable starting point derived from the data if none provided.
    if initial is None:
        w_min = min(tv)
        k_at_min = ks[tv.index(w_min)]
        initial = SVIParams(a=max(w_min, 1e-6), b=0.1, rho=-0.3, m=k_at_min, s=0.1)

    # Unconstrained -> constrained mapping.
    #   b = softplus(pb) >= 0 ; s = softplus(ps) > 0 ; rho = tanh(pr) in (-1,1)
    def softplus(x):
        # Numerically stable.
        return math.log1p(math.exp(-abs(x))) + max(x, 0.0)

    def unpack(p):
        a, pb, pr, m, ps = p
        return SVIParams(a=a, b=softplus(pb), rho=math.tanh(pr), m=m,
                         s=softplus(ps) + 1e-6)

    def pack(sp: SVIParams):
        # Invert softplus/tanh for the initial guess (approximate is fine).
        def inv_softplus(y):
            y = max(y, 1e-9)
            return math.log(math.expm1(y)) if y < 30 else y
        pr = math.atanh(max(min(sp.rho, 0.999), -0.999))
        return [sp.a, inv_softplus(sp.b), pr, sp.m, inv_softplus(max(sp.s - 1e-6, 1e-6))]

    def objective(p):
        sp = unpack(p)
        err = 0.0
        for i in range(n):
            diff = sp.total_variance(ks[i]) - tv[i]
            err += weights[i] * diff * diff
        return err

    # Multi-start: raw SVI has flat valleys where Nelder-Mead can stall on a
    # degenerate (huge-b) fit, so try several seeds and keep the best. Seeds are
    # fixed (no RNG) to stay deterministic across runs.
    w_min = min(tv)
    k_at_min = ks[tv.index(w_min)]
    seeds = [initial]
    for b0 in (0.05, 0.2, 0.5):
        for rho0 in (-0.5, 0.0, 0.3):
            seeds.append(SVIParams(a=max(w_min * 0.5, 1e-6), b=b0, rho=rho0,
                                   m=k_at_min, s=0.2))

    best_p, best_f = None, float("inf")
    for seed in seeds:
        p, f = nelder_mead(objective, pack(seed), step=0.2,
                           max_iter=max_iter, tol=1e-14)
        if f < best_f:
            best_p, best_f = p, f

    params = unpack(best_p)
    rmse = math.sqrt(best_f / wsum)
    return params, rmse
