"""Continuously-monitored lookback options (closed form under BSM).

A lookback option pays off against the extreme (max or min) of the underlying's
path rather than its terminal value:

  * **floating-strike** — the strike is the realized extreme, so the holder
    buys at the lowest (call) or sells at the highest (put) price on the path.
    Priced by Goldman-Sosin-Gatto (1979).
  * **fixed-strike** — an ordinary strike ``K`` applied to the realized extreme
    (call pays ``max(S_max - K, 0)``, put pays ``max(K - S_min, 0)``). Priced
    by Conze-Viswanathan (1991).

Both closed forms take the running extreme observed so far (``s_max`` /
``s_min``); at inception pass the current spot. Formulas assume continuous
monitoring and use cost of carry ``b`` (defaults to ``r``).
"""

import math

from .mathfns import norm_cdf
from .bsm import OptionType, _coerce_type, _validate


def _check(S, t, sigma, b):
    if S <= 0:
        raise ValueError("spot S must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if sigma < 0:
        raise ValueError("sigma must be non-negative")


def floating_strike_lookback(S, t, r, sigma, option_type=OptionType.CALL,
                             s_extreme=None, b=None) -> float:
    """Floating-strike lookback (Goldman-Sosin-Gatto).

    Args:
        s_extreme: running minimum (for a call) or maximum (for a put) observed
            so far. Defaults to the current spot (inception).
        b: cost of carry (defaults to r).

    Call payoff: ``S_T - S_min``. Put payoff: ``S_max - S_T``.
    """
    ot = _coerce_type(option_type)
    _check(S, t, sigma, b)
    if b is None:
        b = r
    if s_extreme is None:
        s_extreme = S
    # The Goldman-Sosin-Gatto formula divides by b; b==0 is a removable
    # singularity, so nudge it by a tiny epsilon (the price is continuous in b).
    if abs(b) < 1e-8:
        b = 1e-8
    carry = math.exp((b - r) * t)
    disc = math.exp(-r * t)

    if t == 0 or sigma == 0:
        if ot is OptionType.CALL:
            return disc * max(S * math.exp(b * t) - s_extreme, 0.0)
        return disc * max(s_extreme - S * math.exp(b * t), 0.0)

    vsqrt = sigma * math.sqrt(t)
    v2 = sigma * sigma
    coef = v2 / (2.0 * b)

    if ot is OptionType.CALL:
        m = s_extreme  # running minimum
        a1 = (math.log(S / m) + (b + 0.5 * v2) * t) / vsqrt
        a2 = a1 - vsqrt
        return (S * carry * norm_cdf(a1)
                - m * disc * norm_cdf(a2)
                + S * disc * coef * (
                    ((S / m) ** (-2.0 * b / v2))
                    * norm_cdf(-a1 + 2.0 * b / sigma * math.sqrt(t))
                    - math.exp(b * t) * norm_cdf(-a1)))
    else:
        M = s_extreme  # running maximum
        b1 = (math.log(S / M) + (b + 0.5 * v2) * t) / vsqrt
        b2 = b1 - vsqrt
        return (M * disc * norm_cdf(-b2)
                - S * carry * norm_cdf(-b1)
                + S * disc * coef * (
                    -((S / M) ** (-2.0 * b / v2))
                    * norm_cdf(b1 - 2.0 * b / sigma * math.sqrt(t))
                    + math.exp(b * t) * norm_cdf(b1)))


def fixed_strike_lookback(S, K, t, r, sigma, option_type=OptionType.CALL,
                          s_extreme=None, b=None) -> float:
    """Fixed-strike lookback (Conze-Viswanathan).

    Call pays ``max(S_max - K, 0)``; put pays ``max(K - S_min, 0)``.

    Args:
        s_extreme: running maximum (call) or minimum (put) so far. Defaults to
            the current spot.
    """
    ot = _coerce_type(option_type)
    _check(S, t, sigma, b)
    if K <= 0:
        raise ValueError("strike K must be positive")
    if b is None:
        b = r
    if s_extreme is None:
        s_extreme = S
    if abs(b) < 1e-8:
        b = 1e-8  # removable singularity in the b-divided coefficient
    carry = math.exp((b - r) * t)
    disc = math.exp(-r * t)

    if t == 0 or sigma == 0:
        if ot is OptionType.CALL:
            return disc * max(max(s_extreme, S * math.exp(b * t)) - K, 0.0)
        return disc * max(K - min(s_extreme, S * math.exp(b * t)), 0.0)

    vsqrt = sigma * math.sqrt(t)
    v2 = sigma * sigma
    coef = v2 / (2.0 * b)

    if ot is OptionType.CALL:
        M = s_extreme  # running max
        if K > M:
            d1 = (math.log(S / K) + (b + 0.5 * v2) * t) / vsqrt
            d2 = d1 - vsqrt
            return (S * carry * norm_cdf(d1) - K * disc * norm_cdf(d2)
                    + S * disc * coef * (
                        -((S / K) ** (-2.0 * b / v2))
                        * norm_cdf(d1 - 2.0 * b / sigma * math.sqrt(t))
                        + math.exp(b * t) * norm_cdf(d1)))
        e1 = (math.log(S / M) + (b + 0.5 * v2) * t) / vsqrt
        e2 = e1 - vsqrt
        return (disc * (M - K)
                + S * carry * norm_cdf(e1) - M * disc * norm_cdf(e2)
                + S * disc * coef * (
                    -((S / M) ** (-2.0 * b / v2))
                    * norm_cdf(e1 - 2.0 * b / sigma * math.sqrt(t))
                    + math.exp(b * t) * norm_cdf(e1)))
    else:
        m = s_extreme  # running min
        if K < m:
            d1 = (math.log(S / K) + (b + 0.5 * v2) * t) / vsqrt
            d2 = d1 - vsqrt
            return (K * disc * norm_cdf(-d2) - S * carry * norm_cdf(-d1)
                    + S * disc * coef * (
                        ((S / K) ** (-2.0 * b / v2))
                        * norm_cdf(-d1 + 2.0 * b / sigma * math.sqrt(t))
                        - math.exp(b * t) * norm_cdf(-d1)))
        f1 = (math.log(S / m) + (b + 0.5 * v2) * t) / vsqrt
        f2 = f1 - vsqrt
        return (disc * (K - m)
                - S * carry * norm_cdf(-f1) + m * disc * norm_cdf(-f2)
                + S * disc * coef * (
                    ((S / m) ** (-2.0 * b / v2))
                    * norm_cdf(-f1 + 2.0 * b / sigma * math.sqrt(t))
                    - math.exp(b * t) * norm_cdf(-f1)))
