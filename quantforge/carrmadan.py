"""Shared Carr-Madan Fourier pricer for exponential-Levy models.

Any model whose log-price is ``ln S_T = ln S + (r - q + omega) t + L_t`` for a
Levy process ``L`` with a known characteristic exponent ``psi`` (so
``E[e^{i u L_t}] = e^{t psi(u)}``) can be priced by the Carr-Madan (1999)
transform: damp the call price by ``e^{alpha k}`` to make it integrable, then
recover it from a single Fourier integral of the risk-neutral characteristic
function. The martingale drift correction is ``omega = -psi(-i)``.

This module factors that machinery out so each Levy model only supplies its
``psi`` (see :mod:`quantforge.cgmy`, :mod:`quantforge.nig`); the Gauss-Legendre
integration nodes are shared with the Heston pricer. Pure standard library.
"""

import cmath
import math

from .bsm import OptionType, _coerce_type
from .heston import _GL_NODES, _GL_WEIGHTS


def carr_madan_call(S, K, t, r, q, psi, alpha=1.5, upper=200.0):
    """Carr-Madan European call price for a Levy model with exponent ``psi``.

    Args:
        psi: callable ``psi(u)`` returning the complex Levy characteristic
            exponent, ``E[e^{i u L_t}] = exp(t * psi(u))``.
        alpha: damping factor (> 0); needs ``E[S_T^{alpha+1}] < inf``, i.e.
            ``psi(-(alpha+1) i)`` finite.
        upper: truncation of the Fourier integral.

    Returns the undiscounted-to-priced European call value.
    """
    x0 = math.log(S)
    lnK = math.log(K)
    disc = math.exp(-r * t)
    omega = -psi(-1j)   # martingale correction: E[S_T] = S e^{(r-q)t}

    def char_logspot(u):
        drift = x0 + (r - q + omega) * t
        return cmath.exp(1j * u * drift + t * psi(u))

    half = 0.5 * upper
    total = 0.0
    for node, w in zip(_GL_NODES, _GL_WEIGHTS):
        nu = half * (node + 1.0)
        if nu <= 0:
            nu = 1e-8
        u = nu - (alpha + 1.0) * 1j
        phi = char_logspot(u)
        denom = alpha * alpha + alpha - nu * nu + 1j * (2.0 * alpha + 1.0) * nu
        rho = disc * phi / denom
        total += w * (cmath.exp(-1j * nu * lnK) * rho).real
    return math.exp(-alpha * lnK) * half * total / math.pi


def levy_price(S, K, t, r, q, psi, option_type=OptionType.CALL,
               alpha=1.5, upper=200.0) -> float:
    """Price a European call/put for a Levy model via Carr-Madan + parity."""
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if t == 0:
        return max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)
    call = carr_madan_call(S, K, t, r, q, psi, alpha=alpha, upper=upper)
    if ot is OptionType.CALL:
        return call
    return call - S * math.exp(-q * t) + K * math.exp(-r * t)
