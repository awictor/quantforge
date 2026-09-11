"""Coupon-bond yield analytics: price, yield-to-maturity, duration, convexity.

These are the standard single-yield fixed-income risk measures, working off an
explicit cashflow schedule ``(time, amount)`` and a continuously-compounded
yield ``y`` unless noted:

    price       = sum_i CF_i e^{-y t_i}
    Macaulay D  = sum_i t_i CF_i e^{-y t_i} / price          (years)
    modified D  = -1/price dPrice/dy = Macaulay D             (continuous comp.)
    convexity   = 1/price d2Price/dy2 = sum_i t_i^2 CF_i e^{-y t_i} / price
    DV01        = -dPrice/dy * 1e-4 = modified D * price * 1e-4

With continuous compounding modified and Macaulay duration coincide. A helper
builds the level-coupon schedule of a vanilla bond, and the yield-to-maturity is
solved from a price by Newton with a bisection fallback. Pure standard library.
"""

import math
from typing import List, Sequence, Tuple


def bond_cashflows(face, coupon_rate, maturity, freq=2) -> List[Tuple[float, float]]:
    """Level-coupon schedule ``[(t, amount), ...]`` for a vanilla bond.

    ``freq`` coupons per year of ``face * coupon_rate / freq`` each, with the
    face repaid alongside the final coupon at ``maturity``.
    """
    if freq < 1:
        raise ValueError("freq must be >= 1")
    if maturity <= 0:
        raise ValueError("maturity must be positive")
    n = int(round(maturity * freq))
    if n < 1:
        raise ValueError("need at least one coupon period")
    cpn = face * coupon_rate / freq
    flows = []
    for i in range(1, n + 1):
        t = i / freq
        amt = cpn + (face if i == n else 0.0)
        flows.append((t, amt))
    return flows


def bond_price_from_yield(cashflows: Sequence[Tuple[float, float]], y) -> float:
    """Present value of the cashflows at continuously-compounded yield ``y``."""
    return sum(cf * math.exp(-y * t) for t, cf in cashflows)


def macaulay_duration(cashflows: Sequence[Tuple[float, float]], y) -> float:
    """Macaulay duration (years): PV-weighted average cashflow time."""
    price = bond_price_from_yield(cashflows, y)
    if price <= 0:
        raise ValueError("bond price must be positive")
    return sum(t * cf * math.exp(-y * t) for t, cf in cashflows) / price


def modified_duration(cashflows: Sequence[Tuple[float, float]], y) -> float:
    """Modified duration ``-1/price dPrice/dy``.

    Under continuous compounding this equals the Macaulay duration.
    """
    return macaulay_duration(cashflows, y)


def convexity(cashflows: Sequence[Tuple[float, float]], y) -> float:
    """Convexity ``1/price d2Price/dy2``: PV-weighted average of squared time."""
    price = bond_price_from_yield(cashflows, y)
    if price <= 0:
        raise ValueError("bond price must be positive")
    return sum(t * t * cf * math.exp(-y * t) for t, cf in cashflows) / price


def dv01(cashflows: Sequence[Tuple[float, float]], y) -> float:
    """Dollar value of a 1bp yield rise (negative: prices fall as yields rise).

    ``DV01 = dPrice/dy * 1e-4 = -modified_duration * price * 1e-4``.
    """
    price = bond_price_from_yield(cashflows, y)
    return -modified_duration(cashflows, y) * price * 1e-4


def price_from_curve(cashflows: Sequence[Tuple[float, float]], curve) -> float:
    """Present value of the cashflows off a discount curve.

    ``curve`` is anything callable as ``curve.df(t)`` (e.g.
    :class:`quantforge.DiscountCurve`) or a plain ``curve(t)`` returning the
    discount factor ``P(0, t)``.
    """
    df = curve.df if hasattr(curve, "df") else curve
    return sum(cf * df(t) for t, cf in cashflows)


def key_rate_durations(cashflows, pillar_times, zero_rates, bump=1e-4):
    """Key-rate (partial) durations of a bond against a zero-rate curve.

    Builds a log-linear :class:`quantforge.DiscountCurve` from the pillar zero
    rates, then bumps each pillar's zero rate up by ``bump`` in turn and measures
    the fractional price change ``-dP/P / bump``. Returns a list of key-rate
    durations aligned with ``pillar_times``; their sum approximates the bond's
    effective duration (a parallel shift is the sum of the pillar bumps).

    Continuously-compounded zero rates; ``DF = e^{-z t}`` at each pillar.
    """
    from .discount_curve import DiscountCurve
    if len(pillar_times) != len(zero_rates):
        raise ValueError("pillar_times and zero_rates must have equal length")
    base_curve = DiscountCurve.from_zero_rates(pillar_times, zero_rates)
    base = price_from_curve(cashflows, base_curve)
    if base <= 0:
        raise ValueError("bond price must be positive")
    krds = []
    for i in range(len(pillar_times)):
        bumped = list(zero_rates)
        bumped[i] += bump
        up = price_from_curve(cashflows,
                              DiscountCurve.from_zero_rates(pillar_times, bumped))
        krds.append(-(up - base) / base / bump)
    return krds


def effective_duration_from_curve(cashflows, pillar_times, zero_rates,
                                  bump=1e-4):
    """Effective duration under a parallel shift of the whole zero curve.

    Shifts every pillar zero rate by +/- ``bump`` and central-differences the
    fractional price change. Equals the sum of the :func:`key_rate_durations` to
    first order.
    """
    from .discount_curve import DiscountCurve
    if len(pillar_times) != len(zero_rates):
        raise ValueError("pillar_times and zero_rates must have equal length")
    up = price_from_curve(
        cashflows, DiscountCurve.from_zero_rates(
            pillar_times, [z + bump for z in zero_rates]))
    dn = price_from_curve(
        cashflows, DiscountCurve.from_zero_rates(
            pillar_times, [z - bump for z in zero_rates]))
    base = price_from_curve(
        cashflows, DiscountCurve.from_zero_rates(pillar_times, zero_rates))
    if base <= 0:
        raise ValueError("bond price must be positive")
    return -(up - dn) / (2.0 * bump) / base


def yield_to_maturity(cashflows: Sequence[Tuple[float, float]], price,
                      tol=1e-10, max_iter=100) -> float:
    """Solve the continuously-compounded yield reproducing ``price``.

    Newton on the price/yield relation (derivative is ``-D * price``) with a
    bracketing bisection fallback, since price is monotone decreasing in yield.
    """
    if price <= 0:
        raise ValueError("price must be positive")
    total = sum(cf for _, cf in cashflows)
    if price > total + 1e-12:
        raise ValueError("price exceeds the sum of undiscounted cashflows")
    lo, hi = -0.5, 5.0  # yields between -50% and 500% (cts. comp.)
    # Price decreases in y: p(lo) is the highest, p(hi) the lowest.
    y = 0.05
    for _ in range(max_iter):
        p = bond_price_from_yield(cashflows, y)
        diff = p - price
        if abs(diff) < tol:
            return y
        deriv = -sum(t * cf * math.exp(-y * t) for t, cf in cashflows)
        if deriv == 0.0:
            break
        step = diff / deriv
        y_new = y - step
        if not (lo < y_new < hi):
            # Fall back to bisection on the maintained bracket.
            if diff > 0:  # model price too high -> raise yield
                lo = y
            else:
                hi = y
            y_new = 0.5 * (lo + hi)
        else:
            if diff > 0:
                lo = y
            else:
                hi = y
        y = y_new
    return y
