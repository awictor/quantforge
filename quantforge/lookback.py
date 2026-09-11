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


# Broadie-Glasserman-Kou (1999) continuity-correction constant,
# beta = -zeta(1/2)/sqrt(2 pi).
_BGK_BETA = 0.5826


def discrete_fixed_strike_lookback(S, K, t, r, sigma, n_fixings,
                                   option_type=OptionType.CALL, b=None) -> float:
    """Discretely-monitored fixed-strike lookback (Broadie-Glasserman-Kou 1999).

    The realized extreme is sampled at ``n_fixings`` equally-spaced dates rather
    than continuously, which lowers a call-on-max and raises a put-on-min versus
    continuous monitoring. Broadie-Glasserman-Kou give an asymptotic continuity
    correction: shift the *spot* fed to the continuous
    :func:`fixed_strike_lookback` by ``exp(-/+ beta sigma sqrt(dt))`` (down for a
    call on the max, up for a put on the min), with ``beta ~ 0.5826`` and
    ``dt = t / n_fixings``. As ``n_fixings -> infinity`` the shift vanishes and
    the price converges to the continuous lookback.

    Accurate to a few tenths of a percent for ``n_fixings`` of ~50 or more; the
    correction is asymptotic, so coarse monitoring (a handful of dates) carries a
    larger error.
    """
    ot = _coerce_type(option_type)
    _check(S, t, sigma, b)
    if K <= 0:
        raise ValueError("strike K must be positive")
    if n_fixings < 1:
        raise ValueError("n_fixings must be >= 1")
    if b is None:
        b = r
    if t == 0 or sigma == 0:
        return fixed_strike_lookback(S, K, t, r, sigma, ot, b=b)
    dt = t / n_fixings
    shift = math.exp(_BGK_BETA * sigma * math.sqrt(dt))
    # Scaling the monitored extreme by `shift` is equivalent to scaling the
    # initial spot in the extreme's distribution: down for a max, up for a min.
    s_eff = S / shift if ot is OptionType.CALL else S * shift
    return fixed_strike_lookback(s_eff, K, t, r, sigma, ot, b=b)


def lookback_greeks(S, t, r, sigma, option_type=OptionType.CALL, b=None,
                    kind="floating", K=None, s_extreme=None):
    """Greeks of a lookback option by central finite differences.

    ``kind`` selects the closed form: ``"floating"``
    (:func:`floating_strike_lookback`) or ``"fixed"``
    (:func:`fixed_strike_lookback`, which needs ``K``). Returns a dict with
    delta, gamma, vega, and theta (calendar, per year). ``s_extreme`` (the
    running min/max) defaults to the current spot.
    """
    ot = _coerce_type(option_type)
    _check(S, t, sigma, b)
    if b is None:
        b = r

    if kind == "floating":
        def px(S_=S, t_=t, sigma_=sigma):
            return floating_strike_lookback(S_, t_, r, sigma_, ot,
                                            s_extreme=s_extreme, b=b)
    elif kind == "fixed":
        if K is None:
            raise ValueError("fixed-strike lookback needs K")

        def px(S_=S, t_=t, sigma_=sigma):
            return fixed_strike_lookback(S_, K, t_, r, sigma_, ot,
                                         s_extreme=s_extreme, b=b)
    else:
        raise ValueError("kind must be 'floating' or 'fixed'")

    base = px()
    hS = 1e-3 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)

    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)

    ht = min(1e-4, 0.5 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)

    return {"price": base, "delta": delta, "gamma": gamma,
            "vega": vega, "theta": theta}


def discrete_fixed_strike_lookback_greeks(S, K, t, r, sigma, n_fixings,
                                          option_type=OptionType.CALL, b=None):
    """Greeks of a discretely-monitored fixed-strike lookback by central finite
    differences of :func:`discrete_fixed_strike_lookback`: ``delta`` (dV/dS),
    ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), ``theta`` (calendar decay). The
    number of monitoring dates ``n_fixings`` is held fixed. Returns a dict with
    ``price`` and those fields.
    """
    ot = _coerce_type(option_type)
    _check(S, t, sigma, b)
    if b is None:
        b = r

    def px(S_=S, t_=t, sigma_=sigma):
        return discrete_fixed_strike_lookback(S_, K, t_, r, sigma_, n_fixings,
                                              ot, b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}
