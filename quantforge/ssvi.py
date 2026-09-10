"""Surface SVI (SSVI): an arbitrage-free whole-surface volatility parametrization.

Gatheral & Jacquier (2014), "Arbitrage-free SVI volatility surfaces", tie every
expiry's smile together through a single skew function ``phi`` and the term
structure of at-the-money total variance ``theta_t = w(0, t)``:

    w(k, theta) = theta/2 * ( 1 + rho * phi(theta) * k
                              + sqrt( (phi(theta) * k + rho)^2 + (1 - rho^2) ) )

Each fixed-``t`` slice is an ordinary SVI smile, but the shared ``(rho, phi)``
makes the surface consistent across maturities and lets simple conditions on
``phi`` rule out both butterfly (within-slice) and calendar (across-slice) static
arbitrage. This module uses the power-law skew

    phi(theta) = eta / ( theta^gamma * (1 + theta)^(1 - gamma) ),   0 < gamma < 1,

fits ``(rho, eta, gamma)`` plus one ``theta`` per expiry to a market vol surface,
and exposes the Gatheral-Jacquier no-arbitrage checks. Pure standard library.
"""

import math
from dataclasses import dataclass, field
from typing import Dict, Sequence, Tuple, List

from .optimize import nelder_mead


def ssvi_phi(theta: float, eta: float, gamma: float) -> float:
    """Power-law SSVI skew function phi(theta)."""
    return eta / (theta ** gamma * (1.0 + theta) ** (1.0 - gamma))


def ssvi_total_variance(k: float, theta: float, rho: float,
                        eta: float, gamma: float) -> float:
    """SSVI total implied variance w(k, theta)."""
    ph = ssvi_phi(theta, eta, gamma)
    x = ph * k
    return 0.5 * theta * (1.0 + rho * x
                          + math.sqrt((x + rho) ** 2 + (1.0 - rho * rho)))


@dataclass(frozen=True)
class SSVIParams:
    """A fitted SSVI surface: global (rho, eta, gamma) and per-expiry theta."""
    rho: float
    eta: float
    gamma: float
    thetas: Dict[float, float] = field(default_factory=dict)

    def total_variance(self, k: float, t: float) -> float:
        return ssvi_total_variance(k, self.thetas[t], self.rho,
                                   self.eta, self.gamma)

    def implied_vol(self, k: float, t: float) -> float:
        w = self.total_variance(k, t)
        return math.sqrt(max(w, 0.0) / t)

    def phi(self, t: float) -> float:
        return ssvi_phi(self.thetas[t], self.eta, self.gamma)


def ssvi_butterfly_free(theta: float, rho: float, eta: float, gamma: float,
                        tol: float = 1e-9) -> bool:
    """Gatheral-Jacquier sufficient condition for a butterfly-arbitrage-free slice.

    A fixed-``theta`` SSVI slice has no butterfly (density-negative) arbitrage if

        theta * phi * (1 + |rho|) <= 4      and
        theta * phi^2 * (1 + |rho|) <= 4 .
    """
    ph = ssvi_phi(theta, eta, gamma)
    c1 = theta * ph * (1.0 + abs(rho))
    c2 = theta * ph * ph * (1.0 + abs(rho))
    return c1 <= 4.0 + tol and c2 <= 4.0 + tol


def ssvi_calendar_free(params: SSVIParams, ks: Sequence[float] = None,
                       tol: float = 1e-9) -> bool:
    """Check the surface has no calendar-spread arbitrage on the fitted expiries.

    Calendar arbitrage is absent when total variance is non-decreasing in
    maturity at every log-moneyness: ``w(k, t_{i+1}) >= w(k, t_i)``. Checked on a
    grid of ``k`` across each adjacent pair of the fitted expiries.
    """
    if ks is None:
        ks = [x * 0.1 for x in range(-15, 16)]  # k in [-1.5, 1.5]
    ts = sorted(params.thetas)
    for t_lo, t_hi in zip(ts, ts[1:]):
        for k in ks:
            w_lo = ssvi_total_variance(k, params.thetas[t_lo], params.rho,
                                       params.eta, params.gamma)
            w_hi = ssvi_total_variance(k, params.thetas[t_hi], params.rho,
                                       params.eta, params.gamma)
            if w_hi < w_lo - tol:
                return False
    return True


def ssvi_is_arbitrage_free(params: SSVIParams, ks: Sequence[float] = None) -> bool:
    """True if every slice is butterfly-free and the surface is calendar-free."""
    for t, theta in params.thetas.items():
        if not ssvi_butterfly_free(theta, params.rho, params.eta, params.gamma):
            return False
    return ssvi_calendar_free(params, ks)


def calibrate_ssvi(
    market: Sequence[Tuple[float, float, float]],
    initial: SSVIParams = None, max_iter: int = 8000,
) -> Tuple[SSVIParams, float]:
    """Fit an SSVI surface to market implied vols.

    ``market`` is a sequence of ``(t, k, iv)`` points (expiry in years,
    log-moneyness, Black-Scholes implied vol). Fits the global ``(rho, eta,
    gamma)`` and one ``theta`` per distinct expiry by minimising the total-
    variance RMSE, using a smooth reparametrization so ``theta > 0``, ``eta > 0``,
    ``gamma in (0, 1)`` and ``rho in (-1, 1)`` and the built-in Nelder-Mead.

    Returns ``(params, rmse)`` where ``rmse`` is the root-mean-square implied-vol
    error across the market points.
    """
    pts = [(float(t), float(k), float(iv)) for t, k, iv in market]
    if len(pts) < 5:
        raise ValueError("need at least five (t, k, iv) points")
    expiries = sorted({t for t, _, _ in pts})
    n_exp = len(expiries)
    idx = {t: i for i, t in enumerate(expiries)}
    # Market total variance targets.
    w_mkt = [(idx[t], k, iv * iv * t) for t, k, iv in pts]

    def softplus(x):
        return math.log1p(math.exp(-abs(x))) + max(x, 0.0)

    def inv_softplus(y):
        y = max(y, 1e-8)
        return math.log(math.expm1(y)) if y < 30 else y

    def sigmoid(x):
        return 1.0 / (1.0 + math.exp(-x))

    def unpack(p):
        rho = math.tanh(p[0])
        eta = softplus(p[1]) + 1e-8
        gamma = min(max(sigmoid(p[2]), 1e-6), 1.0 - 1e-6)
        thetas = [softplus(p[3 + i]) + 1e-10 for i in range(n_exp)]
        return rho, eta, gamma, thetas

    if initial is None:
        # Seed each theta from the ATM-ish market variance at that expiry.
        theta0 = []
        for t in expiries:
            near = min((pt for pt in pts if pt[0] == t), key=lambda pt: abs(pt[1]))
            theta0.append(max(near[2] ** 2 * t, 1e-6))
        eta0, gamma0, rho0 = 1.0, 0.5, -0.3
    else:
        eta0, gamma0, rho0 = initial.eta, initial.gamma, initial.rho
        theta0 = [initial.thetas[t] for t in expiries]

    x0 = [math.atanh(max(min(rho0, 0.999), -0.999)),
          inv_softplus(eta0),
          math.log(gamma0 / (1.0 - gamma0))]
    x0 += [inv_softplus(th) for th in theta0]

    def objective(p):
        rho, eta, gamma, thetas = unpack(p)
        err = 0.0
        for i, k, w in w_mkt:
            model = ssvi_total_variance(k, thetas[i], rho, eta, gamma)
            diff = model - w
            err += diff * diff
        return err

    best_p, _ = nelder_mead(objective, x0, step=0.3, max_iter=max_iter, tol=1e-16)
    rho, eta, gamma, thetas = unpack(best_p)
    params = SSVIParams(rho=rho, eta=eta, gamma=gamma,
                        thetas={t: thetas[idx[t]] for t in expiries})

    # Report RMSE in implied-vol units.
    sse = 0.0
    for t, k, iv in pts:
        w = params.total_variance(k, t)
        model_iv = math.sqrt(max(w, 0.0) / t)
        sse += (model_iv - iv) ** 2
    return params, math.sqrt(sse / len(pts))
