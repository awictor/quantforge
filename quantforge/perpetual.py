"""Perpetual American options (no expiry): exact closed form.

A perpetual American option never expires, so the optimal-exercise problem is
time-homogeneous and has a clean closed form (Merton 1973). The value is a power
function of spot up to a flat exercise boundary, beyond which the holder
exercises for intrinsic value.

For a perpetual call with cost of carry ``b`` and rate ``r`` (both positive, and
``r > b`` so the call is finite):

    h1 = 0.5 - b/sigma^2 + sqrt((b/sigma^2 - 0.5)^2 + 2 r / sigma^2)
    S* = K * h1 / (h1 - 1)                        (exercise boundary)
    C  = (K / (h1 - 1)) * ((h1 - 1)/h1 * S/K)^h1   for S < S*, else S - K

with the symmetric formula for the put. These match the finite-maturity
American price as the maturity grows large, which the tests check against the
binomial tree at long expiry.
"""

import math

from .bsm import OptionType, _coerce_type


def perpetual_american(S, K, r, sigma, option_type=OptionType.CALL, b=None):
    """Price a perpetual American option (Merton 1973), exact closed form.

    Args:
        b: cost of carry (defaults to r). A perpetual call requires ``b < r``
            (some carry cost / dividend) to be finite and worth exercising;
            with ``b >= r`` the call is never exercised early and its value
            tends to the spot.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    if r <= 0:
        raise ValueError("r must be positive for a finite perpetual value")
    if b is None:
        b = r

    v2 = sigma * sigma
    if ot is OptionType.CALL:
        # h1 > 1 root of the characteristic quadratic.
        h1 = 0.5 - b / v2 + math.sqrt((b / v2 - 0.5) ** 2 + 2.0 * r / v2)
        if h1 <= 1.0:
            # b >= r: never optimal to exercise; value approaches the asset.
            return S
        S_star = K * h1 / (h1 - 1.0)
        if S >= S_star:
            return S - K
        return (S_star - K) * (S / S_star) ** h1
    else:
        # h2 < 0 root.
        h2 = 0.5 - b / v2 - math.sqrt((b / v2 - 0.5) ** 2 + 2.0 * r / v2)
        S_star = K * h2 / (h2 - 1.0)
        if S <= S_star:
            return K - S
        return (K - S_star) * (S / S_star) ** h2


def perpetual_exercise_boundary(K, r, sigma, option_type=OptionType.CALL, b=None):
    """The optimal-exercise spot ``S*`` for a perpetual American option."""
    ot = _coerce_type(option_type)
    if b is None:
        b = r
    v2 = sigma * sigma
    if ot is OptionType.CALL:
        h1 = 0.5 - b / v2 + math.sqrt((b / v2 - 0.5) ** 2 + 2.0 * r / v2)
        if h1 <= 1.0:
            return math.inf
        return K * h1 / (h1 - 1.0)
    h2 = 0.5 - b / v2 - math.sqrt((b / v2 - 0.5) ** 2 + 2.0 * r / v2)
    return K * h2 / (h2 - 1.0)
