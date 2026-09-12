"""Capital-budgeting metrics: NPV, IRR, MIRR, payback, profitability index.

Discounted-cashflow project appraisal. Cashflows are given as a list starting
with the (usually negative) time-0 outlay followed by the period cashflows. Pure
standard library.
"""

import math


def npv(rate, cashflows):
    """Net present value at a per-period discount ``rate``.

    ``sum_t CF_t / (1 + rate)^t`` with ``t = 0`` the first entry. Falls as the
    discount rate rises for a conventional project.
    """
    if rate <= -1.0:
        raise ValueError("rate must exceed -100%")
    return sum(cf / (1.0 + rate) ** t for t, cf in enumerate(cashflows))


def irr(cashflows, tol=1e-10, max_iter=200):
    """Internal rate of return: the rate at which :func:`npv` is zero.

    Bisection on ``[-0.999, 10]`` (NPV is monotone decreasing in the rate for a
    conventional outlay-then-inflows project). Requires a sign change in the
    cashflows. Recovers the discount rate that generated a set of flows.
    """
    if not cashflows or all(cf >= 0 for cf in cashflows) or all(cf <= 0 for cf in cashflows):
        raise ValueError("cashflows must contain both signs")
    lo, hi = -0.999, 10.0
    f_lo, f_hi = npv(lo, cashflows), npv(hi, cashflows)
    if f_lo * f_hi > 0:
        raise ValueError("no sign change in NPV over the search range")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        fm = npv(mid, cashflows)
        if abs(fm) < tol:
            return mid
        if f_lo * fm < 0:
            hi = mid
        else:
            lo, f_lo = mid, fm
    return 0.5 * (lo + hi)


def profitability_index(rate, cashflows):
    """Profitability index: PV of inflows over the initial outlay.

    ``PV(CF_1..) / -CF_0`` -- above one for a value-adding project, exactly one at
    the :func:`irr`. Requires a negative time-0 outlay.
    """
    if cashflows[0] >= 0:
        raise ValueError("first cashflow must be a negative outlay")
    pv_inflows = sum(cashflows[t] / (1.0 + rate) ** t
                     for t in range(1, len(cashflows)))
    return pv_inflows / (-cashflows[0])


def payback_period(cashflows):
    """Payback period: fractional periods to recover the initial outlay.

    Accumulates undiscounted cashflows until the running total turns non-negative,
    interpolating within the crossing period. Returns ``inf`` if never recovered.
    """
    if cashflows[0] >= 0:
        raise ValueError("first cashflow must be a negative outlay")
    cum = cashflows[0]
    for t in range(1, len(cashflows)):
        prev = cum
        cum += cashflows[t]
        if cum >= 0:
            # Linear interpolation across period t.
            return (t - 1) + (-prev) / cashflows[t] if cashflows[t] != 0 else float(t)
    return float("inf")


def mirr(cashflows, finance_rate, reinvest_rate):
    """Modified internal rate of return.

    Compounds positive cashflows forward at ``reinvest_rate`` to the terminal date,
    discounts negative cashflows back at ``finance_rate`` to time zero, then
    ``MIRR = (FV_pos / -PV_neg)^{1/n} - 1``. Avoids the multiple-IRR problem and
    uses realistic reinvestment.
    """
    n = len(cashflows) - 1
    if n < 1:
        raise ValueError("need at least two cashflows")
    fv_pos = sum(cf * (1.0 + reinvest_rate) ** (n - t)
                 for t, cf in enumerate(cashflows) if cf > 0)
    pv_neg = sum(cf / (1.0 + finance_rate) ** t
                 for t, cf in enumerate(cashflows) if cf < 0)
    if pv_neg == 0 or fv_pos == 0:
        raise ValueError("need both positive and negative cashflows")
    return (fv_pos / (-pv_neg)) ** (1.0 / n) - 1.0
