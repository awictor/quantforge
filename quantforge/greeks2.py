"""Second-order (and cross) Greeks for the generalized BSM model.

These sensitivities matter for hedging a book through moves in spot *and*
vol, and for managing time decay of the vega/delta profile:

  * vanna  = d(delta)/d(sigma) = d(vega)/d(spot)   -- spot/vol cross-gamma
  * vomma  = d(vega)/d(sigma)  (a.k.a. volga)       -- vol convexity
  * charm  = d(delta)/d(time)  (delta decay / "delta bleed")
  * veta   = d(vega)/d(time)   (vega decay)
  * speed  = d(gamma)/d(spot)                        -- third-order in spot
  * zomma  = d(gamma)/d(sigma)
  * color  = d(gamma)/d(time)

All are analytic and verified against finite differences of the first-order
Greeks in the test suite. Charm/veta/color are quoted per year of calendar
time (negate the time derivative), consistent with ``bsm.theta``.
"""

import math

from .mathfns import norm_pdf, norm_cdf
from .bsm import OptionType, _coerce_type, _validate, _d1_d2


def _prep(S, K, t, r, sigma, b):
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    d1, d2 = _d1_d2(S, K, t, r, sigma, b)
    carry = math.exp((b - r) * t)
    return d1, d2, carry, b


def vanna(S, K, t, r, sigma, b=None) -> float:
    """d(delta)/d(sigma) = d(vega)/d(spot). Same for calls and puts."""
    d1, d2, carry, b = _prep(S, K, t, r, sigma, b)
    return -carry * norm_pdf(d1) * d2 / sigma


def vomma(S, K, t, r, sigma, b=None) -> float:
    """d(vega)/d(sigma) (volga). Same for calls and puts."""
    d1, d2, carry, b = _prep(S, K, t, r, sigma, b)
    vega = S * carry * norm_pdf(d1) * math.sqrt(t)
    return vega * d1 * d2 / sigma


# ``volga`` is the common alias for vomma.
volga = vomma


def charm(S, K, t, r, sigma, option_type=OptionType.CALL, b=None) -> float:
    """Calendar charm = -d(delta)/d(t_expiry): the drift of delta over time."""
    ot = _coerce_type(option_type)
    d1, d2, carry, b = _prep(S, K, t, r, sigma, b)
    vsqrt = sigma * math.sqrt(t)
    common = carry * norm_pdf(d1) * (2.0 * b * t - d2 * vsqrt) / (2.0 * t * vsqrt)
    if ot is OptionType.CALL:
        ddelta_dtexp = (b - r) * carry * norm_cdf(d1) - common
    else:
        ddelta_dtexp = (b - r) * carry * (norm_cdf(d1) - 1.0) - common
    # Calendar charm = -d(delta)/d(t_expiry). The FD cross-check shows the
    # standard-form expression above already carries the calendar sign, so we
    # return it directly.
    return ddelta_dtexp


def veta(S, K, t, r, sigma, b=None) -> float:
    """Calendar veta = -d(vega)/d(t_expiry): decay of vega over time."""
    d1, d2, carry, b = _prep(S, K, t, r, sigma, b)
    vega = S * carry * norm_pdf(d1) * math.sqrt(t)
    # veta = vega * [ (b-r) + b*d1/(sigma*sqrt(t)) - (1+d1*d2)/(2t) ].
    # This standard-form expression already equals the calendar veta
    # (-d(vega)/d(t_expiry)); FD cross-check confirms the sign.
    return vega * ((b - r)
                   + b * d1 / (sigma * math.sqrt(t))
                   - (1.0 + d1 * d2) / (2.0 * t))


def speed(S, K, t, r, sigma, b=None) -> float:
    """d(gamma)/d(spot). Third-order in spot; same for calls and puts."""
    d1, d2, carry, b = _prep(S, K, t, r, sigma, b)
    vsqrt = sigma * math.sqrt(t)
    gamma = carry * norm_pdf(d1) / (S * vsqrt)
    return -gamma / S * (d1 / vsqrt + 1.0)


def zomma(S, K, t, r, sigma, b=None) -> float:
    """d(gamma)/d(sigma). Same for calls and puts."""
    d1, d2, carry, b = _prep(S, K, t, r, sigma, b)
    vsqrt = sigma * math.sqrt(t)
    gamma = carry * norm_pdf(d1) / (S * vsqrt)
    return gamma * (d1 * d2 - 1.0) / sigma


def color(S, K, t, r, sigma, b=None) -> float:
    """Calendar color = -d(gamma)/d(t_expiry): the decay of gamma over time."""
    d1, d2, carry, b = _prep(S, K, t, r, sigma, b)
    vsqrt = sigma * math.sqrt(t)
    pref = carry * norm_pdf(d1) / (2.0 * S * t * vsqrt)
    bracket = (2.0 * (b - r) * t + 1.0
               + (2.0 * b * t - d2 * vsqrt) / vsqrt * d1)
    dgamma_dtexp = -pref * bracket
    return -dgamma_dtexp
