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


def _svi_derivs(p: SVIParams, k):
    """Total variance w(k) and its first two derivatives w'(k), w''(k)."""
    d = k - p.m
    root = math.sqrt(d * d + p.s * p.s)
    w = p.a + p.b * (p.rho * d + root)
    wp = p.b * (p.rho + d / root)
    wpp = p.b * p.s * p.s / (root * root * root)
    return w, wp, wpp


def svi_variance_swap_strike(p: SVIParams, S0, t, r, q=0.0, n_strikes=401,
                             width=8.0):
    """Fair variance-swap strike (annualized *variance*) implied by a raw-SVI slice.

    Replicates the variance swap from the SVI smile: at each strike the Black
    implied vol is ``p.implied_vol(k, t)`` with ``k = ln(K / F)`` the
    log-moneyness on the forward ``F = S0 e^{(r-q)t}``. Feeds the smile to
    :func:`quantforge.variance_swap_from_smile`, so the result is model-
    consistent with the fitted slice. Returns the fair *variance* (square it back
    to vol with ``sqrt``); a flat slice (``b = 0``) returns that flat variance
    ``sigma^2``, and a skewed slice returns a variance above the ATM variance (the
    convexity/skew premium).
    """
    from .varswap import variance_swap_from_smile

    F = S0 * math.exp((r - q) * t)

    def vol_fn(K):
        return p.implied_vol(math.log(K / F), t)

    return variance_swap_from_smile(S0, t, r, vol_fn, q=q, n_strikes=n_strikes,
                                    width=width)


def svi_g(p: SVIParams, k):
    """Gatheral-Jacquier g-function of a raw-SVI slice at log-moneyness ``k``.

    The slice is free of butterfly (static/density) arbitrage iff ``g(k) >= 0``
    for all ``k`` (the implied risk-neutral density is then non-negative). With
    ``w = w(k)``, ``w'`` and ``w''``:

        g(k) = (1 - k w' / (2w))^2 - (w'^2 / 4)(1/w + 1/4) + w''/2.
    """
    w, wp, wpp = _svi_derivs(p, k)
    if w <= 0:
        return float("-inf")
    term1 = (1.0 - k * wp / (2.0 * w)) ** 2
    term2 = (wp * wp / 4.0) * (1.0 / w + 0.25)
    return term1 - term2 + wpp / 2.0


def svi_butterfly_arbitrage(p: SVIParams, ks=None, tol=1e-10):
    """Return the log-moneyness points where the SVI slice has butterfly arb.

    Scans ``ks`` (default a wide grid) and reports those where ``g(k) < -tol``.
    An empty list means the slice is butterfly-arbitrage-free on the grid.
    """
    if ks is None:
        ks = [(-2.0 + 4.0 * i / 400.0) for i in range(401)]
    return [k for k in ks if svi_g(p, k) < -tol]


def svi_is_butterfly_free(p: SVIParams, ks=None) -> bool:
    return not svi_butterfly_arbitrage(p, ks)


def lee_wing_slopes(p: SVIParams):
    """Asymptotic wing slopes of total variance for a raw-SVI slice.

    As ``k -> +/- inf`` the SVI total variance ``w(k)`` is linear with slopes

        right (k -> +inf):  b (1 + rho)
        left  (k -> -inf):  b (1 - rho)

    Lee's moment formula caps the slope of *total variance* at 2 for a valid
    (arbitrage-free-wing) surface, so both slopes must be <= 2. Returns
    ``(left_slope, right_slope)``.
    """
    return (p.b * (1.0 - p.rho), p.b * (1.0 + p.rho))


def lee_bounds_ok(p: SVIParams, tol=1e-9):
    """True if both SVI wing slopes satisfy Lee's ``slope <= 2`` moment bound.

    Equivalent to :meth:`SVIParams.is_arbitrage_free_wings` but exposes the two
    directional slopes explicitly via :func:`lee_wing_slopes`.
    """
    left, right = lee_wing_slopes(p)
    return left <= 2.0 + tol and right <= 2.0 + tol


def svi_local_variance(p: SVIParams, k, dw_dt):
    """Dupire local variance of a single SVI slice, analytic in strike.

    Given the slice ``p`` and the total-variance time derivative ``dw_dt =
    dw/dt`` at log-moneyness ``k`` (supplied by the caller, since one slice
    carries no maturity information), the Gatheral total-variance Dupire formula
    gives

        sigma_loc^2 = dw/dt
            / [ 1 - (k/w) w_k + (1/4)(-1/4 - 1/w + k^2/w^2) w_k^2 + (1/2) w_kk ]

    with ``w``, ``w_k = w'(k)`` and ``w_kk = w''(k)`` taken in closed form from
    the SVI parametrization (no finite differences in strike). Raises if the
    Dupire denominator is non-positive (a butterfly-arbitrage flag: the slice's
    ``svi_g`` is negative there).
    """
    w, wp, wpp = _svi_derivs(p, k)
    if w <= 0.0:
        raise ValueError("total variance must be positive")
    denom = (1.0
             - (k / w) * wp
             + 0.25 * (-0.25 - 1.0 / w + k * k / (w * w)) * wp * wp
             + 0.5 * wpp)
    if denom <= 0.0:
        raise ValueError("non-positive Dupire denominator (butterfly arbitrage)")
    return dw_dt / denom


def svi_surface_local_vol(slices, k, t):
    """Local volatility from a term structure of SVI slices at ``(k, t)``.

    ``slices`` maps expiry ``t_i`` (years) to a fitted :class:`SVIParams`. Total
    variance is interpolated *linearly in t* at fixed ``k`` to supply the
    Dupire ``dw/dt`` (the piecewise-constant slope of the bracketing slices),
    while the strike derivatives come analytically from the slice active at
    ``t``. ``t`` must lie within the fitted expiry range.

    Returns the local volatility ``sqrt(sigma_loc^2)``.
    """
    ts = sorted(slices)
    if not ts:
        raise ValueError("need at least one SVI slice")
    if t < ts[0] - 1e-12 or t > ts[-1] + 1e-12:
        raise ValueError("t outside the fitted expiry range")
    # Bracketing slices for the linear w(k, .) slope.
    if len(ts) == 1:
        raise ValueError("need at least two slices for a time derivative")
    if t <= ts[0]:
        i = 0
    elif t >= ts[-1]:
        i = len(ts) - 2
    else:
        i = max(j for j in range(len(ts) - 1) if ts[j] <= t)
    t0, t1 = ts[i], ts[i + 1]
    p0, p1 = slices[t0], slices[t1]
    w0 = p0.total_variance(k)
    w1 = p1.total_variance(k)
    dw_dt = (w1 - w0) / (t1 - t0)
    # Interpolate the slice-local strike derivatives at t using the same linear
    # weight, so w and its k-derivatives are consistent with the interpolated w.
    lam = (t - t0) / (t1 - t0)
    w_a, wp_a, wpp_a = _svi_derivs(p0, k)
    w_b, wp_b, wpp_b = _svi_derivs(p1, k)
    w = (1 - lam) * w_a + lam * w_b
    wp = (1 - lam) * wp_a + lam * wp_b
    wpp = (1 - lam) * wpp_a + lam * wpp_b
    if w <= 0.0:
        raise ValueError("total variance must be positive")
    denom = (1.0
             - (k / w) * wp
             + 0.25 * (-0.25 - 1.0 / w + k * k / (w * w)) * wp * wp
             + 0.5 * wpp)
    if denom <= 0.0:
        raise ValueError("non-positive Dupire denominator (butterfly arbitrage)")
    return math.sqrt(dw_dt / denom)


def svi_repair_butterfly(p: SVIParams, ks=None, max_iter=200, factor=0.98):
    """Repair a single SVI slice's butterfly arbitrage by shrinking the wings.

    If the slice has ``svi_g(k) < 0`` anywhere (a negative density), it reduces
    the wing angle ``b`` geometrically (which flattens the smile and lifts the
    g-function) until :func:`svi_is_butterfly_free` passes or ``max_iter`` is
    reached. Returns a new :class:`SVIParams`; the ATM level, skew, shift and
    curvature are preserved. If already arbitrage-free the input is returned
    unchanged.
    """
    if svi_is_butterfly_free(p, ks):
        return p
    b = p.b
    for _ in range(max_iter):
        b *= factor
        candidate = SVIParams(a=p.a, b=b, rho=p.rho, m=p.m, s=p.s)
        if svi_is_butterfly_free(candidate, ks):
            return candidate
    return SVIParams(a=p.a, b=b, rho=p.rho, m=p.m, s=p.s)


def calibrate_svi_from_prices(F, t, r, strikes, call_prices, q=0.0,
                              vega_weighted=True, initial=None, max_iter=4000):
    """Calibrate a raw SVI slice directly from market *call prices*.

    Inverts each call to its Black-Scholes implied volatility, converts to total
    variance ``w = sigma^2 t``, and fits raw SVI with :func:`calibrate_svi`.
    Quotes are vega-weighted by default (near-the-money prices carry the most
    volatility information, so weighting by Black vega down-weights the deep
    wings where a price error maps to a large vol error).

    Args:
        F: forward. strikes, call_prices: matching market quotes at expiry ``t``.
        r: discount rate (the calls are priced on the forward, carry ``b = r``
            relative to spot ``S = F e^{-rt}``... here calls are taken on the
            forward directly with discounting ``e^{-rt}``).
        vega_weighted: weight each quote by its Black vega if True.

    Returns ``(params, iv_rmse, price_rmse)``.
    """
    from .implied import implied_volatility
    from .bsm import call_price as bs_call, vega as bs_vega, OptionType

    strikes = [float(k) for k in strikes]
    call_prices = [float(c) for c in call_prices]
    n = len(strikes)
    if n != len(call_prices) or n < 5:
        raise ValueError("need at least five matching (strike, call) quotes")

    S = F * math.exp(-r * t)   # spot consistent with the given forward
    ks, tv, weights = [], [], []
    ivs = []
    for K, c in zip(strikes, call_prices):
        try:
            iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=r)
        except ValueError:
            continue
        ks.append(math.log(K / F))
        tv.append(iv * iv * t)
        ivs.append((K, iv))
        w = bs_vega(S, K, t, r, iv, b=r) if vega_weighted else 1.0
        weights.append(max(w, 1e-8))
    if len(ks) < 5:
        raise ValueError("fewer than five prices could be inverted to vols")

    params, iv_rmse = calibrate_svi(ks, tv, weights=weights, initial=initial,
                                    max_iter=max_iter)

    # Price RMSE: reprice each strike at the fitted SVI vol.
    sse = 0.0
    for K, c in zip(strikes, call_prices):
        k = math.log(K / F)
        w = params.total_variance(k)
        if w <= 0:
            continue
        iv = math.sqrt(w / t)
        model_c = bs_call(S, K, t, r, iv, b=r)
        sse += (model_c - c) ** 2
    price_rmse = math.sqrt(sse / n)
    return params, iv_rmse, price_rmse
