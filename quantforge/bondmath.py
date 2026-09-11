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


def dated_bond_cashflows(start, maturity_years, face, coupon_rate, freq=2,
                         convention="30/360", end_of_month=False,
                         business_day="unadjusted"):
    """Coupon-bond cashflows on a real calendar with day-count accruals.

    Generates the coupon schedule from ``start`` (a ``(y, m, d)`` date) with
    :func:`quantforge.generate_schedule`, computes each period's accrual factor
    with :func:`quantforge.year_fraction` under ``convention``, and pays the
    day-count-weighted coupon ``face * coupon_rate * tau_i`` each period plus the
    face at maturity. Returns ``[(pay_date, year_fraction, amount), ...]``.

    Unlike :func:`bond_cashflows` (which assumes uniform ``1/freq`` periods),
    this reflects the actual day counts, month-end roll, and business-day
    adjustment -- the difference that matters for act/360 and stub periods.
    """
    from .schedule import generate_schedule
    from .daycount import year_fraction
    if coupon_rate < 0:
        raise ValueError("coupon_rate must be non-negative")
    dates = generate_schedule(start, maturity_years, 12 // freq
                              if 12 % freq == 0 else max(1, round(12 / freq)),
                              end_of_month=end_of_month,
                              convention=business_day)
    flows = []
    prev = start
    n = len(dates)
    for i, pay in enumerate(dates):
        tau = year_fraction(prev, pay, convention)
        cpn = face * coupon_rate * tau
        amt = cpn + (face if i == n - 1 else 0.0)
        flows.append((pay, tau, amt))
        prev = pay
    return flows


def accrued_interest(face, coupon_rate, freq, fraction_elapsed) -> float:
    """Accrued interest since the last coupon, straight-line within the period.

    ``fraction_elapsed`` in ``[0, 1]`` is the share of the current coupon period
    that has passed at settlement. The accrual is
    ``face * coupon_rate / freq * fraction_elapsed`` -- the linear (actual/
    actual-in-period) convention. The buyer pays this on top of the quoted
    clean price.
    """
    if not (0.0 <= fraction_elapsed <= 1.0):
        raise ValueError("fraction_elapsed must be in [0, 1]")
    if freq < 1:
        raise ValueError("freq must be >= 1")
    return face * coupon_rate / freq * fraction_elapsed


def dirty_price(cashflows: Sequence[Tuple[float, float]], y) -> float:
    """Dirty (invoice) price: the full present value of the remaining cashflows.

    Alias of :func:`bond_price_from_yield` -- the cash amount actually paid at
    settlement, before subtracting accrued interest to get the clean quote.
    """
    return bond_price_from_yield(cashflows, y)


def clean_price(cashflows: Sequence[Tuple[float, float]], y, face,
                coupon_rate, freq, fraction_elapsed) -> float:
    """Clean (quoted) price: dirty price minus accrued interest.

    ``clean = dirty - accrued``. At a coupon date (``fraction_elapsed = 0``) the
    clean and dirty prices coincide; mid-period the clean price strips out the
    accrued coupon so the quote does not saw-tooth across coupon dates.
    """
    dirty = bond_price_from_yield(cashflows, y)
    accrued = accrued_interest(face, coupon_rate, freq, fraction_elapsed)
    return dirty - accrued


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
