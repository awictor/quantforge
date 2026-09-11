"""Barone-Adesi-Whaley (1987) quadratic approximation for American options.

BAW splits the American price into the European value plus an early-exercise
premium approximated by the dominant solution of the (time-independent)
quadratic accompanying ODE. The premium is a power of the spot anchored at a
critical price ``S*`` (the exercise boundary at t=0), found by a one-dimensional
Newton solve of the value-matching / smooth-pasting condition. The result is a
fast closed-form approximation that is accurate to a few cents for short and
medium maturities and is the classic benchmark the Kim integral equation and the
binomial tree refine.

Pure standard library; carry ``b = r - q`` (dividend yield ``q``).
"""

import math

from .bsm import OptionType, _coerce_type, _validate, call_price, put_price
from .mathfns import norm_cdf, norm_pdf


def _d1(S, K, t, r, b, sigma):
    return (math.log(S / K) + (b + 0.5 * sigma * sigma) * t) / (sigma * math.sqrt(t))


def _baw_call(S, K, t, r, sigma, b):
    """BAW American call price (b < r, i.e. positive dividend, else European)."""
    if b >= r:
        return call_price(S, K, t, r, sigma, b=b)   # never exercise early

    M = 2.0 * r / (sigma * sigma)
    N = 2.0 * b / (sigma * sigma)
    K_fac = 1.0 - math.exp(-r * t)
    q2 = (-(N - 1.0) + math.sqrt((N - 1.0) ** 2 + 4.0 * M / K_fac)) / 2.0

    Sx = _critical_call(S, K, t, r, sigma, b, q2)
    if S >= Sx:
        return S - K
    d1 = _d1(Sx, K, t, r, b, sigma)
    A2 = (Sx / q2) * (1.0 - math.exp((b - r) * t) * norm_cdf(d1))
    return call_price(S, K, t, r, sigma, b=b) + A2 * (S / Sx) ** q2


def _critical_call(S, K, t, r, sigma, b, q2):
    """Newton solve for the American-call critical spot S* (value matching)."""
    Sx = _seed_critical_call(S, K, t, r, sigma, b)
    for _ in range(100):
        d1 = _d1(Sx, K, t, r, b, sigma)
        eur = call_price(Sx, K, t, r, sigma, b=b)
        carry = math.exp((b - r) * t)
        lhs = eur + (1.0 - carry * norm_cdf(d1)) * Sx / q2
        rhs = Sx - K
        bi = (carry * norm_cdf(d1) * (1.0 - 1.0 / q2)
              + (1.0 - carry * norm_pdf(d1) / (sigma * math.sqrt(t))) / q2)
        f = lhs - rhs
        if abs(f) < 1e-8:
            break
        Sx -= f / (bi - 1.0)
        if Sx <= 0:
            Sx = 1e-6 * K
    return Sx


def _seed_critical_call(S, K, t, r, sigma, b):
    """Starting guess for the American-call critical spot (BAW seed)."""
    M = 2.0 * r / (sigma * sigma)
    N = 2.0 * b / (sigma * sigma)
    q2u = (-(N - 1.0) + math.sqrt((N - 1.0) ** 2 + 4.0 * M)) / 2.0
    Su = K / (1.0 - 1.0 / q2u)
    h2 = -(b * t + 2.0 * sigma * math.sqrt(t)) * K / (Su - K)
    return K + (Su - K) * (1.0 - math.exp(h2))


def _baw_put(S, K, t, r, sigma, b):
    """BAW American put price."""
    M = 2.0 * r / (sigma * sigma)
    N = 2.0 * b / (sigma * sigma)
    K_fac = 1.0 - math.exp(-r * t)
    q1 = (-(N - 1.0) - math.sqrt((N - 1.0) ** 2 + 4.0 * M / K_fac)) / 2.0

    Sx = _critical_put(S, K, t, r, sigma, b, q1)
    if S <= Sx:
        return K - S
    d1 = _d1(Sx, K, t, r, b, sigma)
    A1 = -(Sx / q1) * (1.0 - math.exp((b - r) * t) * norm_cdf(-d1))
    return put_price(S, K, t, r, sigma, b=b) + A1 * (S / Sx) ** q1


def _critical_put(S, K, t, r, sigma, b, q1):
    """Newton solve for the American-put critical spot S* (value matching)."""
    Sx = _seed_critical_put(S, K, t, r, sigma, b)
    for _ in range(100):
        d1 = _d1(Sx, K, t, r, b, sigma)
        eur = put_price(Sx, K, t, r, sigma, b=b)
        carry = math.exp((b - r) * t)
        lhs = eur - (1.0 - carry * norm_cdf(-d1)) * Sx / q1
        rhs = K - Sx
        bi = (-carry * norm_cdf(-d1) * (1.0 - 1.0 / q1)
              - (1.0 + carry * norm_pdf(-d1) / (sigma * math.sqrt(t))) / q1)
        f = lhs - rhs
        if abs(f) < 1e-8:
            break
        Sx -= f / (bi + 1.0)
        if Sx <= 0:
            Sx = 1e-6 * K
    return Sx


def _seed_critical_put(S, K, t, r, sigma, b):
    """Starting guess for the American-put critical spot (BAW seed)."""
    M = 2.0 * r / (sigma * sigma)
    N = 2.0 * b / (sigma * sigma)
    q1u = (-(N - 1.0) - math.sqrt((N - 1.0) ** 2 + 4.0 * M)) / 2.0
    Su = K / (1.0 - 1.0 / q1u)
    h1 = (b * t - 2.0 * sigma * math.sqrt(t)) * K / (K - Su)
    return Su + (K - Su) * math.exp(h1)


def baw_american(S, K, t, r, sigma, option_type=OptionType.CALL, b=None) -> float:
    """American option price by the Barone-Adesi-Whaley quadratic approximation.

    ``b`` is the cost of carry (defaults to ``r``); dividend yield ``q`` enters
    as ``b = r - q``. A no-dividend American call (``b = r``) returns the
    European value. Falls back to intrinsic below/above the critical spot.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if t == 0:
        return max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)
    if ot is OptionType.CALL:
        return _baw_call(S, K, t, r, sigma, b)
    return _baw_put(S, K, t, r, sigma, b)


def baw_critical_spot(K, t, r, sigma, option_type=OptionType.CALL, b=None):
    """Barone-Adesi-Whaley early-exercise boundary S* at inception.

    The spot at which immediate exercise becomes optimal: a call is exercised
    for ``S >= S*`` and a put for ``S <= S*``. Returns ``None`` when early
    exercise is never optimal (an American call with ``b >= r`` equals its
    European value, so there is no finite boundary).
    """
    ot = _coerce_type(option_type)
    _validate(1.0, K, t, sigma)
    if b is None:
        b = r
    v2 = sigma * sigma
    M = 2.0 * r / v2
    N = 2.0 * b / v2
    K_fac = 1.0 - math.exp(-r * t)
    if ot is OptionType.CALL:
        if b >= r:
            return None  # never exercise a no-dividend American call early
        q2 = (-(N - 1.0) + math.sqrt((N - 1.0) ** 2 + 4.0 * M / K_fac)) / 2.0
        return _critical_call(K, K, t, r, sigma, b, q2)
    q1 = (-(N - 1.0) - math.sqrt((N - 1.0) ** 2 + 4.0 * M / K_fac)) / 2.0
    return _critical_put(K, K, t, r, sigma, b, q1)


def baw_american_greeks(S, K, t, r, sigma, option_type=OptionType.CALL, b=None):
    """Greeks of a BAW American option by central finite differences of
    :func:`baw_american`: ``delta`` (dV/dS), ``gamma`` (d2V/dS2), ``vega``
    (dV/dsigma), ``theta`` (calendar decay). Returns a dict with ``price`` and
    those fields.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r

    def px(S_=S, t_=t, sigma_=sigma):
        return baw_american(S_, K, t_, r, sigma_, ot, b=b)

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
