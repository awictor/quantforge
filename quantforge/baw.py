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

    # Solve for the critical spot S* by Newton on the value-matching equation.
    Sx = _seed_critical_call(S, K, t, r, sigma, b)
    for _ in range(100):
        d1 = _d1(Sx, K, t, r, b, sigma)
        eur = call_price(Sx, K, t, r, sigma, b=b)
        carry = math.exp((b - r) * t)
        lhs = eur + (1.0 - carry * norm_cdf(d1)) * Sx / q2
        rhs = Sx - K
        # Derivative of (lhs - rhs) w.r.t. Sx.
        bi = (carry * norm_cdf(d1) * (1.0 - 1.0 / q2)
              + (1.0 - carry * norm_pdf(d1) / (sigma * math.sqrt(t))) / q2)
        f = lhs - rhs
        if abs(f) < 1e-8:
            break
        Sx -= f / (bi - 1.0)
        if Sx <= 0:
            Sx = 1e-6 * K

    if S >= Sx:
        return S - K
    d1 = _d1(Sx, K, t, r, b, sigma)
    A2 = (Sx / q2) * (1.0 - math.exp((b - r) * t) * norm_cdf(d1))
    return call_price(S, K, t, r, sigma, b=b) + A2 * (S / Sx) ** q2


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

    if S <= Sx:
        return K - S
    d1 = _d1(Sx, K, t, r, b, sigma)
    A1 = -(Sx / q1) * (1.0 - math.exp((b - r) * t) * norm_cdf(-d1))
    return put_price(S, K, t, r, sigma, b=b) + A1 * (S / Sx) ** q1


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
