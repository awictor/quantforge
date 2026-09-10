"""Local-stochastic-volatility (LSV) leverage-function calibration.

An LSV model multiplies a stochastic-variance process by a deterministic
*leverage function* ``L(S, t)`` chosen so the model reprices the market's local
volatility exactly:

    dS/S = L(S, t) sqrt(V_t) dW,   dV = kappa(theta - V) dt + xi sqrt(V) dZ.

Gyongy's theorem gives the calibration condition -- the LSV and the pure
local-vol model share vanilla prices iff

    L(K, t)^2 = sigma_Dupire(K, t)^2 / E[ V_t | S_t = K ].

The conditional expectation ``E[V_t | S_t = K]`` has no closed form, so it is
estimated by simulating the (Heston) variance process and binning ``V_t`` by the
spot level -- the particle method of Guyon & Henry-Labordere (2012). This module
calibrates ``L`` on a time-then-space grid by that fixed point and exposes the
resulting leverage surface. Pure standard library.
"""

import math
import random

from .heston_mc import _norm_ppf


def calibrate_leverage(S0, r, local_vol_fn, kappa, theta, xi, rho, v0,
                       expiries, k_grid, q=0.0, n_paths=20000, seed=None,
                       sub_steps=20):
    """Calibrate the LSV leverage surface by the particle method.

    Args:
        local_vol_fn: target Dupire local vol ``sigma(K, t)``.
        kappa, theta, xi, rho, v0: the backbone Heston variance parameters.
        expiries: increasing calibration times (the leverage is piecewise
            constant in time between them; ``t=0`` implied at the front).
        k_grid: log-moneyness bin edges (relative to the forward) for the
            conditional-expectation binning.

    Returns ``leverage`` -- a dict ``{t: {k_center: L}}`` -- and the callable
    ``lev_fn(spot, t)`` that interpolates it, suitable for an LSV Monte Carlo.
    """
    rng = random.Random(seed)
    n = n_paths
    # Simulate spot and variance jointly, applying the leverage calibrated so
    # far. Leverage at t_0 (short end) is local_vol / sqrt(v0) at the money.
    logS = [math.log(S0)] * n
    V = [v0] * n
    z_perp = math.sqrt(max(1.0 - rho * rho, 0.0))

    leverage = {}
    centers = [0.5 * (k_grid[i] + k_grid[i + 1]) for i in range(len(k_grid) - 1)]

    prev_t = 0.0
    for t in expiries:
        dt = t - prev_t
        F = S0 * math.exp((r - q) * t)
        # Bin current V by log-moneyness to estimate E[V | S].
        sums = [0.0] * len(centers)
        counts = [0] * len(centers)
        for p in range(n):
            k = logS[p] - math.log(F)
            b = _bin_index(k, k_grid)
            if b is not None:
                sums[b] += V[p]
                counts[b] += 1
        lev_t = {}
        for i, kc in enumerate(centers):
            cond_V = sums[i] / counts[i] if counts[i] > 0 else v0
            K = F * math.exp(kc)
            loc = local_vol_fn(K, t if t > 0 else expiries[0])
            lev_t[kc] = loc / math.sqrt(cond_V) if cond_V > 1e-12 else 1.0
        leverage[t] = lev_t

        # Evolve to the next expiry with several Euler sub-steps (a single step
        # over a wide pillar gap is far too coarse and biases the calibration).
        h = dt / sub_steps
        sh = math.sqrt(h)
        for _ in range(sub_steps):
            for p in range(n):
                k = logS[p] - math.log(F)
                L = _interp_center(centers, lev_t, k)
                vol = L * math.sqrt(max(V[p], 0.0))
                z1 = rng.gauss(0.0, 1.0)
                z2 = rng.gauss(0.0, 1.0)
                dW = sh * z1
                dZ = sh * (rho * z1 + z_perp * z2)
                logS[p] += (r - q - 0.5 * vol * vol) * h + vol * dW
                V[p] = abs(V[p] + kappa * (theta - V[p]) * h
                           + xi * math.sqrt(max(V[p], 0.0)) * dZ)
        prev_t = t

    ts = sorted(leverage)

    def lev_fn(spot, tt):
        # Nearest calibrated expiry, interpolate in log-moneyness.
        t = min(ts, key=lambda x: abs(x - tt))
        F = S0 * math.exp((r - q) * t)
        k = math.log(spot / F)
        return _interp_center(centers, leverage[t], k)

    return leverage, lev_fn


def _bin_index(k, edges):
    if k < edges[0] or k >= edges[-1]:
        return None
    lo, hi = 0, len(edges) - 1
    while lo < hi - 1:
        mid = (lo + hi) // 2
        if edges[mid] <= k:
            lo = mid
        else:
            hi = mid
    return lo


def _interp_center(centers, lev_t, k):
    if k <= centers[0]:
        return lev_t[centers[0]]
    if k >= centers[-1]:
        return lev_t[centers[-1]]
    for i in range(len(centers) - 1):
        if centers[i] <= k < centers[i + 1]:
            w = (k - centers[i]) / (centers[i + 1] - centers[i])
            return (1 - w) * lev_t[centers[i]] + w * lev_t[centers[i + 1]]
    return lev_t[centers[-1]]
