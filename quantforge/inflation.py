"""Inflation-linked bond and breakeven-inflation analytics.

Real (TIPS-style) instruments scale their cashflows by the ratio of the
reference price index at settlement to the index at issue. This module handles
the index ratio, the inflation-adjusted principal, the Fisher relation between
nominal and real rates, and the breakeven inflation implied by a nominal/real
yield pair. Pure standard library.
"""

import math


def index_ratio(index_settle, index_base) -> float:
    """Index ratio ``CPI_settle / CPI_base`` used to inflate the principal.

    Above 1 when the price index has risen since issue. Both indices must be
    positive.
    """
    if index_base <= 0 or index_settle <= 0:
        raise ValueError("index levels must be positive")
    return index_settle / index_base


def inflation_adjusted_principal(face, index_settle, index_base) -> float:
    """Inflation-adjusted principal ``face * index_ratio`` (the linker notional)."""
    return face * index_ratio(index_settle, index_base)


def fisher_real_rate(nominal, inflation) -> float:
    """Exact Fisher real rate ``(1 + nominal)/(1 + inflation) - 1``.

    The rate that, compounded with inflation, reproduces the nominal rate. For
    small rates it is approximately ``nominal - inflation``.
    """
    return (1.0 + nominal) / (1.0 + inflation) - 1.0


def fisher_nominal_rate(real, inflation) -> float:
    """Exact Fisher nominal rate ``(1 + real)(1 + inflation) - 1``."""
    return (1.0 + real) * (1.0 + inflation) - 1.0


def breakeven_inflation(nominal_yield, real_yield) -> float:
    """Breakeven inflation implied by a nominal and a real yield (Fisher).

    ``(1 + nominal)/(1 + real) - 1`` -- the inflation rate at which a nominal and
    an inflation-linked bond of the same maturity have equal return. The market's
    inflation expectation (plus risk premium).
    """
    if real_yield <= -1.0:
        raise ValueError("real_yield must exceed -100%")
    return (1.0 + nominal_yield) / (1.0 + real_yield) - 1.0


def real_from_breakeven(nominal_yield, breakeven) -> float:
    """Real yield implied by a nominal yield and a breakeven inflation rate.

    ``(1 + nominal)/(1 + breakeven) - 1`` -- inverse of
    :func:`breakeven_inflation`.
    """
    if breakeven <= -1.0:
        raise ValueError("breakeven must exceed -100%")
    return (1.0 + nominal_yield) / (1.0 + breakeven) - 1.0


def deflation_floored_redemption(face, index_settle, index_base) -> float:
    """TIPS-style redemption with the deflation floor: principal never below par.

    Real (TIPS) principal redeems at ``face * max(index_ratio, 1)`` -- the index
    ratio inflates the principal in inflation, but a cumulative deflation over the
    bond's life cannot pull the redemption below the original face. Equals
    :func:`inflation_adjusted_principal` whenever the index has risen since issue
    (ratio >= 1), and is floored to ``face`` otherwise.
    """
    return face * max(index_ratio(index_settle, index_base), 1.0)


def deflation_floor_value(face, index_settle, index_base) -> float:
    """Intrinsic value of the deflation floor: floored redemption minus unfloored.

    ``face * (max(ratio, 1) - ratio)`` -- zero when the index has risen (the floor
    is out of the money), positive under net deflation. The realized payoff of the
    embedded floor option, ignoring optionality/time value.
    """
    ratio = index_ratio(index_settle, index_base)
    return face * (max(ratio, 1.0) - ratio)


def yoy_inflation_rate(index_prev, index_curr) -> float:
    """Year-on-year inflation ``index_curr / index_prev - 1`` between two fixings."""
    if index_prev <= 0 or index_curr <= 0:
        raise ValueError("index levels must be positive")
    return index_curr / index_prev - 1.0


def zc_inflation_swap_rate(index_start, index_end, years) -> float:
    """Fair annualized rate of a zero-coupon inflation swap.

    A ZC inflation swap exchanges ``(1 + k)^T - 1`` (fixed) for the realized index
    growth ``I_T / I_0 - 1`` (float) at maturity. The par fixed rate that zeroes
    the swap is the annualized index growth ``(I_T / I_0)^(1/T) - 1``, so that
    ``(1 + k)^T * I_0 == I_T`` (the compounding identity tested against).
    """
    if index_start <= 0 or index_end <= 0:
        raise ValueError("index levels must be positive")
    if years <= 0:
        raise ValueError("years must be positive")
    return (index_end / index_start) ** (1.0 / years) - 1.0


def zc_inflation_swap_value(notional, fixed_rate, index_start, index_end,
                            years, discount_factor=1.0) -> float:
    """Value of the inflation leg minus the fixed leg of a ZC inflation swap.

    Inflation-leg receiver's value: ``notional * (I_T/I_0 - (1+k)^T)`` at maturity,
    discounted by ``discount_factor``. Zero at the par :func:`zc_inflation_swap_rate`.
    """
    if index_start <= 0 or index_end <= 0:
        raise ValueError("index levels must be positive")
    if years <= 0:
        raise ValueError("years must be positive")
    realized = index_end / index_start
    fixed = (1.0 + fixed_rate) ** years
    return notional * (realized - fixed) * discount_factor


def inflation_curve_from_zc_swaps(index_base, tenors, zc_rates):
    """Projected index levels implied by a strip of zero-coupon swap rates.

    A ZC inflation swap of maturity ``T`` with fair rate ``k_T`` pins the forward
    index to ``I_0 * (1 + k_T)^T`` (the :func:`zc_inflation_swap_rate` identity).
    Given quotes ``(tenors, zc_rates)`` this returns the matching forward index
    levels ``[I_0 (1 + k_T)^T for T in tenors]`` -- the market-implied inflation
    curve, expressed as projected index fixings. By construction reinverting each
    level through :func:`zc_inflation_swap_rate` recovers the input ``zc_rates``.
    """
    if index_base <= 0:
        raise ValueError("index_base must be positive")
    if len(tenors) != len(zc_rates):
        raise ValueError("tenors and zc_rates must have equal length")
    levels = []
    for T, k in zip(tenors, zc_rates):
        if T <= 0:
            raise ValueError("tenors must be positive")
        levels.append(index_base * (1.0 + k) ** T)
    return levels


def forward_inflation_rate(index_start, index_end, t_start, t_end):
    """Annualized forward inflation between two curve horizons.

    ``(I_end / I_start)^(1/(t_end - t_start)) - 1`` -- the constant annual rate
    linking two projected index levels. Chains with the near-leg rate so that
    ``(1 + spot)^t_start (1 + fwd)^(t_end - t_start) = (1 + spot_end)^t_end`` (the
    no-arbitrage forward/spot relation tested against).
    """
    if index_start <= 0 or index_end <= 0:
        raise ValueError("index levels must be positive")
    dt = t_end - t_start
    if dt <= 0:
        raise ValueError("t_end must exceed t_start")
    return (index_end / index_start) ** (1.0 / dt) - 1.0


def yoy_swap_value(notional, fixed_rate, index_levels, discount_factors,
                   index_prev):
    """Value of a year-on-year inflation swap off a projected index curve.

    Each period ``i`` exchanges the realized year-on-year inflation
    ``I_i / I_{i-1} - 1`` (float, received) for ``fixed_rate`` (paid), on
    ``notional``, discounted by ``discount_factors[i]``. ``index_prev`` is the
    fixing one period before the first ``index_levels`` entry (the base for the
    first YoY ratio). Returns the inflation-receiver's value

        notional * sum_i (I_i/I_{i-1} - 1 - fixed_rate) * DF_i.

    Unlike the single-payment ZC swap this pays the annual inflation each period.
    """
    if len(index_levels) != len(discount_factors):
        raise ValueError("index_levels and discount_factors must have equal length")
    if index_prev <= 0:
        raise ValueError("index_prev must be positive")
    prev = index_prev
    value = 0.0
    for I, df in zip(index_levels, discount_factors):
        if I <= 0:
            raise ValueError("index levels must be positive")
        yoy = I / prev - 1.0
        value += (yoy - fixed_rate) * df
        prev = I
    return notional * value


def reference_cpi(cpi_month_start, cpi_next_month, day, days_in_month):
    """Daily reference index by linear interpolation between two monthly fixings.

    Inflation-linked bonds accrue off a *reference index* that interpolates
    linearly within the month between the anchor CPI for the first of the month
    and the first of the next month (the standard linker daily-indexation rule,
    applied to the lagged CPIs). For settlement on the ``day``-th of a month with
    ``days_in_month`` days:

        ref = cpi_month_start + (day - 1)/days_in_month
                  * (cpi_next_month - cpi_month_start)

    Equals ``cpi_month_start`` on the 1st and approaches ``cpi_next_month`` at
    month end.
    """
    if cpi_month_start <= 0 or cpi_next_month <= 0:
        raise ValueError("CPI levels must be positive")
    if days_in_month <= 0:
        raise ValueError("days_in_month must be positive")
    if not (1 <= day <= days_in_month):
        raise ValueError("day must be in [1, days_in_month]")
    frac = (day - 1) / days_in_month
    return cpi_month_start + frac * (cpi_next_month - cpi_month_start)


def index_ratio_interpolated(cpi_month_start, cpi_next_month, day,
                             days_in_month, cpi_base):
    """Index ratio using the daily-interpolated :func:`reference_cpi`.

    ``reference_cpi(...) / cpi_base`` -- the ratio a linker actually applies to
    its principal on a mid-month settlement, versus the month-boundary
    :func:`index_ratio` which ignores intra-month accrual.
    """
    ref = reference_cpi(cpi_month_start, cpi_next_month, day, days_in_month)
    return index_ratio(ref, cpi_base)


def linker_price(real_cashflows, real_yield, index_settle, index_base) -> float:
    """Dirty price of an inflation-linked bond off real cashflows.

    ``real_cashflows`` is ``[(t, real_amount), ...]`` in constant (issue-date)
    money. Each flow is discounted at the continuously-compounded ``real_yield``
    and then the whole bond is inflated by the settlement index ratio:

        price = (index_settle / index_base) * sum_i real_amount_i e^{-r t_i}

    Because the index ratio multiplies every flow, the price is degree-one
    homogeneous in it -- stripping the ratio recovers a standard real-yield bond
    price (the invariant tested against :mod:`quantforge.bondmath`).
    """
    ratio = index_ratio(index_settle, index_base)
    pv = sum(amt * math.exp(-real_yield * t) for t, amt in real_cashflows)
    return ratio * pv


def linker_real_duration(real_cashflows, real_yield) -> float:
    """Modified duration of a linker w.r.t. its real yield (years).

    ``-1/P dP/dr``. Since the settlement index ratio multiplies the whole price it
    cancels in the fractional sensitivity, so the real duration is the PV-weighted
    average cashflow time of the *real* cashflows -- identical to the standard
    :func:`quantforge.modified_duration` on those flows, independent of the index
    level.
    """
    price = sum(amt * math.exp(-real_yield * t) for t, amt in real_cashflows)
    if price <= 0:
        raise ValueError("real cashflow price must be positive")
    return sum(t * amt * math.exp(-real_yield * t) for t, amt in real_cashflows) / price


def linker_real_convexity(real_cashflows, real_yield) -> float:
    """Convexity of a linker w.r.t. its real yield ``1/P d2P/dr2``.

    PV-weighted average of squared cashflow time on the real cashflows; the index
    ratio cancels, matching :func:`quantforge.convexity`.
    """
    price = sum(amt * math.exp(-real_yield * t) for t, amt in real_cashflows)
    if price <= 0:
        raise ValueError("real cashflow price must be positive")
    return sum(t * t * amt * math.exp(-real_yield * t)
               for t, amt in real_cashflows) / price


def linker_real_dv01(real_cashflows, real_yield, index_settle, index_base) -> float:
    """Dollar value of a 1bp real-yield rise for a linker (negative).

    ``dP/dr * 1e-4 = -duration * price * 1e-4`` on the inflated (dirty) price, so
    unlike the fractional duration this DOES scale with the index ratio.
    """
    price = linker_price(real_cashflows, real_yield, index_settle, index_base)
    return -linker_real_duration(real_cashflows, real_yield) * price * 1e-4


def linker_real_yield(real_cashflows, price, index_settle, index_base,
                      tol=1e-10, max_iter=100) -> float:
    """Continuously-compounded real yield reproducing a linker ``price``.

    Deflates the quoted price by the index ratio and solves the standard real-
    cashflow bond yield by bisection (price is monotone decreasing in the yield).
    Inverse of :func:`linker_price`.
    """
    if price <= 0:
        raise ValueError("price must be positive")
    ratio = index_ratio(index_settle, index_base)
    target = price / ratio  # real (deflated) price

    def real_px(r):
        return sum(amt * math.exp(-r * t) for t, amt in real_cashflows)

    lo, hi = -0.5, 5.0
    p_lo, p_hi = real_px(lo), real_px(hi)
    if not (min(p_lo, p_hi) - 1e-9 <= target <= max(p_lo, p_hi) + 1e-9):
        raise ValueError("price outside the achievable yield range")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        pm = real_px(mid)
        if abs(pm - target) < tol:
            return mid
        if pm > target:   # price too high -> raise yield
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
