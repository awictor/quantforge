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


def weighted_average_life(cashflows, balance):
    """Weighted-average life (years) from projected principal cashflows.

    ``WAL = sum_m (month/12) * total_principal_m / balance``. Uses the
    ``total_principal`` column (index 4) of :func:`mbs_cashflows`. Falls as
    prepayment speeds up (principal returns sooner).
    """
    if balance <= 0:
        raise ValueError("balance must be positive")
    return sum((row[0] / 12.0) * row[4] for row in cashflows) / balance
