"""Displaced-diffusion (shifted lognormal) option pricing.

Rubinstein's displaced-diffusion model assumes ``S + shift`` is lognormal
rather than ``S`` itself. The shift lets the underlying reach values as low as
``-shift`` (so negative rates/prices are allowed) and produces a monotone skew:
a positive shift flattens the smile toward normal-like (Bachelier) behavior,
while ``shift = 0`` recovers Black-Scholes.

Pricing is exact: it is a Black-Scholes price on the displaced variables
``S' = S + shift`` and ``K' = K + shift`` with a rescaled volatility
``sigma' = sigma * S / (S + shift)`` chosen so the at-the-money instantaneous
volatility matches ``sigma`` (a common calibration convention).
"""

import math

from .bsm import price as bsm_price, OptionType, _coerce_type


def displaced_diffusion_price(S, K, t, r, sigma, shift=0.0,
                              option_type=OptionType.CALL, b=None) -> float:
    """Price a European option under the displaced-diffusion model.

    Args:
        shift: the displacement added to spot and strike. ``shift = 0`` is
            Black-Scholes; larger positive shifts push toward normal-model
            behavior and allow the underlying to fall below zero (down to
            ``-shift``).
        b: cost of carry (defaults to r).

    The payoff is unchanged (``max(S_T - K, 0)`` etc.); only the diffusion is
    displaced, so the price equals a BSM price on ``S + shift`` / ``K + shift``.
    """
    ot = _coerce_type(option_type)
    if t < 0:
        raise ValueError("t must be non-negative")
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    if b is None:
        b = r
    # Displaced diffusion permits negative spot/strike (down to -shift); only
    # the shifted variables must be positive.
    if S + shift <= 0 or K + shift <= 0:
        raise ValueError("S + shift and K + shift must be positive")

    S_shift = S + shift
    K_shift = K + shift
    # Rescale vol so the ATM instantaneous vol of the displaced process matches
    # sigma at the current spot: sigma_displaced = sigma * S / (S + shift).
    sigma_d = sigma * S / S_shift

    # The displaced forward must equal the true forward plus the (carried)
    # shift, so use a carry that reproduces F' = S*e^{bt} + shift on S_shift.
    fwd = S * math.exp(b * t)
    fwd_shift = fwd + shift
    b_d = math.log(fwd_shift / S_shift) / t if t > 0 else b

    return bsm_price(S_shift, K_shift, t, r, sigma_d, ot, b=b_d)
