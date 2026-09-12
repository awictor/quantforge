"""Bond carry and roll-down return decomposition.

The expected return of holding a bond over a horizon, assuming the yield curve is
unchanged, splits into *carry* (coupon income plus the pull-to-par / financing
term) and *roll-down* (the price gain from the bond aging down a
upward-sloping curve to a lower yield). This module computes each piece and their
total from a zero-rate curve callable. Pure standard library on top of
:mod:`quantforge.bondmath`.
"""

import math


def _price(cashflows, y):
    return sum(cf * math.exp(-y * t) for t, cf in cashflows)


def _yield_from_curve(curve, t):
    """Zero rate at ``t`` from a curve object (``.zero_rate``) or callable."""
    if hasattr(curve, "zero_rate"):
        return curve.zero_rate(t)
    return curve(t)


def carry_return(coupon_rate, yield_now, horizon, financing_rate=0.0):
    """Carry over a horizon: coupon income plus financing, per unit face.

    ``(coupon_rate - financing_rate) * horizon`` -- the running yield earned net
    of the cost of funding the position, holding prices fixed. Positive when the
    coupon exceeds the financing rate.
    """
    return (coupon_rate - financing_rate) * horizon


def rolldown_return(cashflows, curve, horizon):
    """Roll-down return: the price gain purely from the yield rolling down the curve.

    Isolates the yield-change effect from the time-value growth. Values the
    surviving cashflows (those maturing after the horizon), each at its *shortened*
    maturity ``t - horizon``, under two curves: the rolled yield ``y(t - horizon)``
    versus the unchanged-maturity yield ``y(t)``. The fractional difference is the
    roll-down -- zero on a flat curve (the yield does not change as the bond rolls)
    and positive on an upward-sloping curve (the bond rolls to a lower yield).
    """
    if horizon <= 0:
        raise ValueError("horizon must be positive")
    price_base = 0.0
    price_rolled = 0.0
    for t, cf in cashflows:
        rem = t - horizon
        if rem <= 0:
            continue   # cashflow paid within the horizon (income, not roll)
        price_base += cf * math.exp(-_yield_from_curve(curve, t) * rem)
        price_rolled += cf * math.exp(-_yield_from_curve(curve, rem) * rem)
    if price_base <= 0:
        raise ValueError("base price must be positive")
    return price_rolled / price_base - 1.0


def total_carry_rolldown(cashflows, curve, coupon_rate, horizon,
                         financing_rate=0.0):
    """Total expected return = carry + roll-down over the horizon.

    Sums :func:`carry_return` and :func:`rolldown_return`. The expected holding-
    period return if the curve is unchanged; the standard relative-value carry-
    and-roll number.
    """
    yield_now = _yield_from_curve(curve, cashflows[-1][0])
    carry = carry_return(coupon_rate, yield_now, horizon, financing_rate)
    roll = rolldown_return(cashflows, curve, horizon)
    return carry + roll
