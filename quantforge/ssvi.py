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


def _ssvi_k_derivs(k, theta, rho, eta, gamma):
    """w and its k-derivatives at fixed theta: returns (w, dw_dk, d2w_dk2)."""
    ph = ssvi_phi(theta, eta, gamma)
    x = ph * k
    R = math.sqrt((x + rho) ** 2 + (1.0 - rho * rho))
    w = 0.5 * theta * (1.0 + rho * x + R)
    dw_dk = 0.5 * theta * ph * (rho + (x + rho) / R)
    d2w_dk2 = 0.5 * theta * ph * ph * (1.0 - rho * rho) / (R ** 3)
    return w, dw_dk, d2w_dk2


def _ssvi_dw_dtheta(k, theta, rho, eta, gamma):
    """Partial derivative of total variance w.r.t. theta (holding t's k fixed)."""
    ph = ssvi_phi(theta, eta, gamma)
    dph_dtheta = ph * (-gamma / theta - (1.0 - gamma) / (1.0 + theta))
    x = ph * k
    R = math.sqrt((x + rho) ** 2 + (1.0 - rho * rho))
    # d w / d theta at fixed k, accounting for phi(theta).
    base = 0.5 * (1.0 + rho * x + R)
    dx = dph_dtheta * k
    chain = 0.5 * theta * (rho * dx + (x + rho) / R * dx)
    return base + chain


def ssvi_local_variance(k, t, theta, dtheta_dt, rho, eta, gamma):
    """Dupire local variance of an SSVI surface, fully analytic (Gatheral).

    Given the ATM total variance ``theta = theta(t)`` and its time derivative
    ``dtheta_dt = theta'(t)`` at maturity ``t``, the local variance at
    log-moneyness ``k`` is

        sigma_loc^2 = (dw/dt)
            / [ 1 - (k/w) w_k + (1/4)(-1/4 - 1/w + k^2/w^2) w_k^2 + (1/2) w_kk ]

    with all ``w`` derivatives taken in closed form from the SSVI parametrization
    (no finite differences). ``dw/dt = (dw/dtheta) * theta'(t)``.
    """
    w, w_k, w_kk = _ssvi_k_derivs(k, theta, rho, eta, gamma)
    if w <= 0.0:
        raise ValueError("total variance must be positive")
    dw_dt = _ssvi_dw_dtheta(k, theta, rho, eta, gamma) * dtheta_dt
    denom = (1.0
             - (k / w) * w_k
             + 0.25 * (-0.25 - 1.0 / w + k * k / (w * w)) * w_k * w_k
             + 0.5 * w_kk)
    if denom <= 0.0:
        raise ValueError("non-positive Dupire denominator (butterfly arbitrage)")
    return dw_dt / denom


def ssvi_local_vol_from_params(params: SSVIParams, k, t):
    """Local volatility of a fitted SSVI surface at ``(k, t)``.

    Builds ``theta(t)`` and ``theta'(t)`` by linear interpolation of the fitted
    per-expiry ATM total variances (piecewise-linear in ``t``), then applies the
    analytic :func:`ssvi_local_variance`. ``t`` must lie within the fitted expiry
    range.
    """
    ts = sorted(params.thetas)
    if t <= 0.0:
        raise ValueError("t must be positive")
    if t > ts[-1] + 1e-12:
        raise ValueError("t beyond the longest fitted expiry")
    # Locate the bracketing pillars for a linear theta(t) and its slope.
    if t <= ts[0]:
        i = 0
    elif t >= ts[-1]:
        i = len(ts) - 2
    else:
        i = max(j for j in range(len(ts) - 1) if ts[j] <= t)
    t0, t1 = ts[i], ts[i + 1]
    th0, th1 = params.thetas[t0], params.thetas[t1]
    slope = (th1 - th0) / (t1 - t0)
    theta = th0 + slope * (t - t0)
    return math.sqrt(ssvi_local_variance(k, t, theta, slope,
                                         params.rho, params.eta, params.gamma))


def ssvi_local_vol_fn(params: SSVIParams, S0, r, q=0.0):
    """Build a ``(spot, tau) -> local vol`` callable from a fitted SSVI surface.

    Converts the running spot and elapsed time into the SSVI log-moneyness
    ``k = log(spot / F_tau)`` on the forward ``F_tau = S0 e^{(r - q) tau}`` and
    returns the analytic Dupire local vol there. Suitable as the ``local_vol_fn``
    argument of :func:`quantforge.local_vol_mc`, which simulates the surface.

    Below the shortest fitted expiry the ATM total variance is linearly
    extrapolated toward the origin (``theta -> 0`` as ``tau -> 0``, matching
    SSVI's small-time behaviour) rather than frozen, which is what a Monte Carlo
    path integrating ``tau`` from 0 needs; the long end is clamped to the last
    fitted maturity.
    """
    ts = sorted(params.thetas)
    t_max = ts[-1]

    def lv(spot, tau):
        tt = min(max(tau, 1e-6), t_max)
        F = S0 * math.exp((r - q) * tt)
        k = math.log(spot / F)
        try:
            return ssvi_local_vol_from_params(params, k, tt)
        except ValueError:
            # Degenerate denominator far in the wing: fall back to the local ATM
            # vol so the simulation stays finite.
            return ssvi_local_vol_from_params(params, 0.0, tt)

    return lv


def ssvi_reprice_mc(params: SSVIParams, S0, K, t, r, q=0.0,
                    option_type=None, n_steps=100, n_paths=60_000,
                    antithetic=True, seed=None):
    """Monte Carlo a vanilla under the SSVI local-vol surface it calibrates to.

    Closes the calibrate -> local-vol -> reprice loop: simulates the analytic
    SSVI Dupire surface (via :func:`quantforge.local_vol_mc`) and returns the
    option's :class:`~quantforge.MCResult`. A correct local-vol construction
    reprices the SSVI *implied* smile, so this Monte Carlo price should match the
    closed-form Black-Scholes price at the SSVI implied vol for that strike.
    """
    from .montecarlo import local_vol_mc
    from .bsm import OptionType as _OT

    ot = _OT.CALL if option_type is None else option_type
    lv = ssvi_local_vol_fn(params, S0, r, q)
    return local_vol_mc(S0, K, t, r, lv, option_type=ot, q=q,
                        n_steps=n_steps, n_paths=n_paths,
                        antithetic=antithetic, seed=seed)


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
