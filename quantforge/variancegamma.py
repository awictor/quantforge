"""Variance-Gamma (Madan-Carr-Chang 1998) option pricing.

The Variance-Gamma model replaces the Brownian diffusion with Brownian motion
evaluated at a random, gamma-distributed business time. It is a pure-jump
process with three parameters:

    sigma : volatility of the Brownian component
    nu    : variance rate of the gamma time change (controls kurtosis)
    theta : drift of the Brownian component (controls skew; theta < 0 -> the
            equity-style left skew)

As ``nu -> 0`` the model collapses to Black-Scholes. Pricing uses the VG
characteristic function integrated with the same probability decomposition and
Gauss-Legendre quadrature as :mod:`quantforge.heston` (no SciPy).
"""

import cmath
import math

from .bsm import OptionType, _coerce_type
from .heston import _GL_NODES, _GL_WEIGHTS


def _vg_char(u, S, K, t, r, q, sigma, nu, theta):
    """VG characteristic function of log(S_t), evaluated at complex ``u``."""
    x = math.log(S)
    # Martingale drift correction so E[S_t] = S e^{(r-q)t}.
    omega = math.log(1.0 - theta * nu - 0.5 * sigma * sigma * nu) / nu
    mu = (r - q + omega)
    # phi(u) = exp(i u (x + mu t)) * (1 - i theta nu u + 0.5 sigma^2 nu u^2)^{-t/nu}
    drift = cmath.exp(1j * u * (x + mu * t))
    base = 1.0 - 1j * theta * nu * u + 0.5 * sigma * sigma * nu * u * u
    return drift * base ** (-t / nu)


def _probability(S, K, t, r, q, sigma, nu, theta, j, upper=200.0):
    """P_j via the Heston-style characteristic-function integral for VG."""
    x = math.log(S)
    lnK = math.log(K)
    half = 0.5 * upper
    total = 0.0
    for node, w in zip(_GL_NODES, _GL_WEIGHTS):
        phi = half * (node + 1.0)
        if phi <= 0:
            phi = 1e-8
        if j == 1:
            # Share measure: divide the cf by the forward E[S_t].
            cf = (_vg_char(phi - 1j, S, K, t, r, q, sigma, nu, theta)
                  / _vg_char(-1j, S, K, t, r, q, sigma, nu, theta))
        else:
            cf = _vg_char(phi, S, K, t, r, q, sigma, nu, theta)
        integrand = (cmath.exp(-1j * phi * lnK) * cf / (1j * phi)).real
        total += w * integrand
    return 0.5 + half * total / math.pi


def variance_gamma_price(S, K, t, r, sigma, nu, theta,
                         option_type=OptionType.CALL, q=0.0, upper=200.0):
    """Price a European option under the Variance-Gamma model.

    Args:
        sigma: Brownian volatility. nu: gamma-time variance rate (> 0).
        theta: Brownian drift (skew; negative for an equity left skew).
        q: continuous dividend yield.

    Puts follow from put-call parity. As ``nu -> 0`` the price approaches the
    Black-Scholes value.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if nu <= 0 or sigma <= 0:
        raise ValueError("nu and sigma must be positive")
    # The gamma time-change requires 1 - theta nu - 0.5 sigma^2 nu > 0.
    if 1.0 - theta * nu - 0.5 * sigma * sigma * nu <= 0:
        raise ValueError("parameters violate the VG martingale condition")

    if t == 0:
        intrinsic = max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)
        return intrinsic

    P1 = _probability(S, K, t, r, q, sigma, nu, theta, 1, upper)
    P2 = _probability(S, K, t, r, q, sigma, nu, theta, 2, upper)
    call = S * math.exp(-q * t) * P1 - K * math.exp(-r * t) * P2
    if ot is OptionType.CALL:
        return call
    return call - S * math.exp(-q * t) + K * math.exp(-r * t)
