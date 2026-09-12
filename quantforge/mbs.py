"""Mortgage-backed security cashflows and prepayment conventions.

A fixed-rate mortgage pays a level monthly payment that fully amortizes the
balance over its term. Pools of mortgages prepay, quoted as a conditional
prepayment rate (CPR, annualized) or its monthly equivalent (single monthly
mortality, SMM), often on the PSA ramp (0.2% CPR in month 1 rising by 0.2% per
month to 6% CPR at month 30, flat thereafter). This module builds the level
payment, the amortization schedule, the SMM/CPR conversions, the PSA schedule,
the projected cashflows with prepayment, and the weighted-average life. Pure
standard library.
"""

import math


def monthly_payment(balance, annual_rate, term_months):
    """Level fully-amortizing monthly payment for a fixed-rate mortgage.

    ``P = B * i / (1 - (1 + i)^{-n})`` with monthly rate ``i = annual_rate / 12``.
    At zero rate this is the straight-line ``balance / term_months``.
    """
    if balance <= 0:
        raise ValueError("balance must be positive")
    if term_months < 1:
        raise ValueError("term_months must be a positive integer")
    i = annual_rate / 12.0
    if i == 0.0:
        return balance / term_months
    return balance * i / (1.0 - (1.0 + i) ** (-term_months))


def cpr_to_smm(cpr):
    """Single monthly mortality from an annual CPR: ``1 - (1 - CPR)^{1/12}``."""
    if not (0.0 <= cpr < 1.0):
        raise ValueError("cpr must be in [0, 1)")
    return 1.0 - (1.0 - cpr) ** (1.0 / 12.0)


def smm_to_cpr(smm):
    """Annual CPR from a single monthly mortality: ``1 - (1 - SMM)^{12}``.

    Inverse of :func:`cpr_to_smm`.
    """
    if not (0.0 <= smm < 1.0):
        raise ValueError("smm must be in [0, 1)")
    return 1.0 - (1.0 - smm) ** 12


def psa_cpr(month, psa=100.0):
    """CPR on the PSA ramp at a given loan age, for ``psa`` percent of the model.

    Standard 100 PSA: ``CPR = 0.06 * min(month, 30) / 30`` (0.2%/month ramp to 6%
    at month 30, flat after). Scaled by ``psa / 100`` for other speeds.
    """
    if month < 1:
        raise ValueError("month must be a positive integer")
    base = 0.06 * min(month, 30) / 30.0
    return base * psa / 100.0


def amortization_schedule(balance, annual_rate, term_months):
    """Level-payment amortization schedule with no prepayment.

    Returns ``[(month, interest, principal, ending_balance), ...]``. The balance
    amortizes to (floating-point) zero at the final month.
    """
    pay = monthly_payment(balance, annual_rate, term_months)
    i = annual_rate / 12.0
    bal = balance
    rows = []
    for m in range(1, term_months + 1):
        interest = bal * i
        principal = pay - interest
        if m == term_months:
            principal = bal  # clear any residual rounding at maturity
        bal -= principal
        rows.append((m, interest, principal, max(bal, 0.0)))
    return rows


def mbs_cashflows(balance, annual_rate, term_months, smm=0.0):
    """Projected MBS cashflows with a constant SMM prepayment.

    Each month pays scheduled interest and principal on the surviving balance,
    plus a prepayment of ``smm`` times the balance remaining after the scheduled
    principal. Returns ``[(month, interest, scheduled_principal, prepayment,
    total_principal, ending_balance), ...]``. With ``smm = 0`` the total principal
    matches :func:`amortization_schedule`.
    """
    if not (0.0 <= smm < 1.0):
        raise ValueError("smm must be in [0, 1)")
    i = annual_rate / 12.0
    bal = balance
    rows = []
    for m in range(1, term_months + 1):
        if bal <= 0.0:
            break
        # Payment recomputed on the surviving balance and remaining term.
        pay = monthly_payment(bal, annual_rate, term_months - m + 1)
        interest = bal * i
        sched_principal = min(pay - interest, bal)
        after_sched = bal - sched_principal
        prepay = smm * after_sched
        total_principal = sched_principal + prepay
        bal = after_sched - prepay
        rows.append((m, interest, sched_principal, prepay, total_principal,
                     max(bal, 0.0)))
    return rows


def mbs_cashflows_psa(balance, annual_rate, term_months, psa=100.0):
    """Projected MBS cashflows on the PSA prepayment ramp (age-varying SMM).

    Like :func:`mbs_cashflows` but the monthly prepayment uses the age-dependent
    :func:`psa_cpr` converted to SMM at each month, rather than a constant SMM.
    Returns the same ``[(month, interest, scheduled_principal, prepayment,
    total_principal, ending_balance), ...]`` rows. At ``psa = 0`` it reduces to the
    no-prepayment schedule.
    """
    if psa < 0:
        raise ValueError("psa must be non-negative")
    i = annual_rate / 12.0
    bal = balance
    rows = []
    for m in range(1, term_months + 1):
        if bal <= 0.0:
            break
        smm = cpr_to_smm(psa_cpr(m, psa))
        pay = monthly_payment(bal, annual_rate, term_months - m + 1)
        interest = bal * i
        sched_principal = min(pay - interest, bal)
        after_sched = bal - sched_principal
        prepay = smm * after_sched
        total_principal = sched_principal + prepay
        bal = after_sched - prepay
        rows.append((m, interest, sched_principal, prepay, total_principal,
                     max(bal, 0.0)))
    return rows


def mbs_price(cashflows, annual_yield):
    """Present value of projected MBS cashflows at a monthly-compounded yield.

    ``cashflows`` are :func:`mbs_cashflows` rows; each month's cash is
    ``interest + total_principal`` discounted by ``(1 + y/12)^{-month}``. Monotone
    decreasing in ``annual_yield``.
    """
    i = annual_yield / 12.0
    pv = 0.0
    for row in cashflows:
        m, interest, _sched, _prepay, total_principal, _bal = row
        cash = interest + total_principal
        pv += cash / (1.0 + i) ** m
    return pv


def mbs_yield(cashflows, price, tol=1e-10, max_iter=100):
    """Monthly-compounded annual yield reproducing an MBS ``price``.

    Bisection on :func:`mbs_price` (monotone decreasing in yield). Inverse of
    :func:`mbs_price`.
    """
    if price <= 0:
        raise ValueError("price must be positive")
    lo, hi = -0.5, 5.0
    p_lo, p_hi = mbs_price(cashflows, lo), mbs_price(cashflows, hi)
    if not (min(p_lo, p_hi) - 1e-6 <= price <= max(p_lo, p_hi) + 1e-6):
        raise ValueError("price outside the achievable yield range")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        pm = mbs_price(cashflows, mid)
        if abs(pm - price) < tol:
            return mid
        if pm > price:   # price too high -> raise yield
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def weighted_average_life(cashflows, balance):
    """Weighted-average life (years) from projected principal cashflows.

    ``WAL = sum_m (month/12) * total_principal_m / balance``. Uses the
    ``total_principal`` column (index 4) of :func:`mbs_cashflows`. Falls as
    prepayment speeds up (principal returns sooner).
    """
    if balance <= 0:
        raise ValueError("balance must be positive")
    return sum((row[0] / 12.0) * row[4] for row in cashflows) / balance
