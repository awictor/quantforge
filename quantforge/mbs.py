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


def mbs_price_with_spread(cashflows, zero_rates, spread):
    """Present value discounting each cashflow at its zero rate plus a spread.

    ``zero_rates[i]`` is the monthly-compounded annualized zero rate for the cash
    at ``cashflows[i]``'s month; every flow is discounted at ``zero_rate + spread``
    (a parallel add-on, the static/Z-spread convention). Reduces to
    :func:`mbs_price` at a flat curve.
    """
    if len(zero_rates) != len(cashflows):
        raise ValueError("zero_rates and cashflows must align")
    pv = 0.0
    for row, z in zip(cashflows, zero_rates):
        m, interest, _s, _p, total_principal, _b = row
        cash = interest + total_principal
        i = (z + spread) / 12.0
        pv += cash / (1.0 + i) ** m
    return pv


def mbs_zspread(cashflows, zero_rates, price, tol=1e-12, max_iter=100):
    """Static (Z-) spread over the zero curve reproducing an MBS ``price``.

    Bisection on the constant spread added to every zero rate (price is monotone
    decreasing in the spread). Inverse of :func:`mbs_price_with_spread`.
    """
    if price <= 0:
        raise ValueError("price must be positive")
    lo, hi = -0.5, 5.0
    p_lo, p_hi = (mbs_price_with_spread(cashflows, zero_rates, lo),
                  mbs_price_with_spread(cashflows, zero_rates, hi))
    if not (min(p_lo, p_hi) - 1e-6 <= price <= max(p_lo, p_hi) + 1e-6):
        raise ValueError("price outside the achievable spread range")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        pm = mbs_price_with_spread(cashflows, zero_rates, mid)
        if abs(pm - price) < tol:
            return mid
        if pm > price:   # price too high -> raise spread
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def mbs_effective_duration(cashflows, annual_yield, bump=1e-4):
    """Effective duration of an MBS from a parallel yield bump (central difference).

    ``-(P(y+h) - P(y-h)) / (2 h P(y))`` on :func:`mbs_price`. Assumes the cashflows
    are held fixed (a static-duration measure; true option-adjusted duration would
    re-project prepayment at each bumped yield).
    """
    up = mbs_price(cashflows, annual_yield + bump)
    dn = mbs_price(cashflows, annual_yield - bump)
    base = mbs_price(cashflows, annual_yield)
    if base <= 0:
        raise ValueError("base price must be positive")
    return -(up - dn) / (2.0 * bump * base)


def mbs_effective_convexity(cashflows, annual_yield, bump=1e-4):
    """Effective convexity of an MBS from a parallel yield bump.

    ``(P(y+h) - 2 P(y) + P(y-h)) / (h^2 P(y))`` on :func:`mbs_price` (static
    cashflows).
    """
    up = mbs_price(cashflows, annual_yield + bump)
    dn = mbs_price(cashflows, annual_yield - bump)
    base = mbs_price(cashflows, annual_yield)
    if base <= 0:
        raise ValueError("base price must be positive")
    return (up - 2.0 * base + dn) / (bump * bump * base)


def pac_schedule(balance, annual_rate, term_months, psa_low, psa_high):
    """Planned-amortization-class principal schedule from a PSA collar.

    A PAC bond promises the principal that is available under *both* ends of a PSA
    speed band: at each month the scheduled PAC principal is the minimum of the
    total principal produced at ``psa_low`` and at ``psa_high``
    (:func:`mbs_cashflows_psa`). Returns ``[(month, pac_principal), ...]``. Because
    it is a lower envelope, the PAC schedule is stable for any prepayment speed
    inside the collar -- the support (companion) tranche absorbs the difference.
    """
    if psa_low > psa_high:
        raise ValueError("psa_low must not exceed psa_high")
    low = mbs_cashflows_psa(balance, annual_rate, term_months, psa_low)
    high = mbs_cashflows_psa(balance, annual_rate, term_months, psa_high)
    n = max(len(low), len(high))
    sched = []
    for m in range(n):
        p_low = low[m][4] if m < len(low) else 0.0
        p_high = high[m][4] if m < len(high) else 0.0
        sched.append((m + 1, min(p_low, p_high)))
    return sched


def pac_support_split(cashflows, pac_sched):
    """Allocate pool principal between a PAC band and its support tranche.

    At each month the PAC receives its scheduled principal (from
    :func:`pac_schedule`), capped by what the pool actually produces and by the
    PAC's remaining balance; the support tranche receives the remainder. Any PAC
    shortfall in a slow month is made up from later principal before the support is
    paid. Returns ``(pac_rows, support_rows)`` as ``[(month, principal), ...]``.
    The two principal streams sum to the pool principal each month.
    """
    pac_by_month = {m: p for m, p in pac_sched}
    pac_arrears = 0.0
    pac_rows = []
    support_rows = []
    for row in cashflows:
        m = row[0]
        avail = row[4]
        target = pac_by_month.get(m, 0.0) + pac_arrears
        pac_pay = min(target, avail)
        pac_arrears = target - pac_pay  # unmet PAC principal carried forward
        support_pay = avail - pac_pay
        pac_rows.append((m, pac_pay))
        support_rows.append((m, support_pay))
    return pac_rows, support_rows


def sequential_cmo(cashflows, tranche_sizes):
    """Split MBS principal across sequential (plain-vanilla) CMO tranches.

    Principal from ``cashflows`` (:func:`mbs_cashflows` rows) is paid to tranches
    strictly in order: tranche 0 receives all principal until retired, then
    tranche 1, and so on. ``tranche_sizes`` are the initial tranche balances (must
    sum to the pool's total principal). Returns a list, one per tranche, of
    ``[(month, principal, ending_balance), ...]`` rows. Each tranche's principal
    sums to its size; earlier tranches retire first (shorter WAL).
    """
    total_principal = sum(row[4] for row in cashflows)
    if abs(sum(tranche_sizes) - total_principal) > 1e-2:
        raise ValueError("tranche_sizes must sum to the pool's total principal")
    if any(sz < 0 for sz in tranche_sizes):
        raise ValueError("tranche sizes must be non-negative")
    balances = list(tranche_sizes)
    out = [[] for _ in tranche_sizes]
    active = 0
    for row in cashflows:
        principal = row[4]
        m = row[0]
        # Waterfall the month's principal down the tranche stack.
        for k in range(len(balances)):
            if principal <= 0.0:
                break
            if balances[k] <= 0.0:
                continue
            pay = min(principal, balances[k])
            balances[k] -= pay
            principal -= pay
            out[k].append((m, pay, max(balances[k], 0.0)))
    return out


def tranche_wal(tranche_rows, tranche_size):
    """Weighted-average life (years) of a single CMO tranche.

    ``sum_m (month/12) * principal_m / tranche_size`` over the tranche's principal
    rows from :func:`sequential_cmo`.
    """
    if tranche_size <= 0:
        raise ValueError("tranche_size must be positive")
    return sum((row[0] / 12.0) * row[1] for row in tranche_rows) / tranche_size


def weighted_average_life(cashflows, balance):
    """Weighted-average life (years) from projected principal cashflows.

    ``WAL = sum_m (month/12) * total_principal_m / balance``. Uses the
    ``total_principal`` column (index 4) of :func:`mbs_cashflows`. Falls as
    prepayment speeds up (principal returns sooner).
    """
    if balance <= 0:
        raise ValueError("balance must be positive")
    return sum((row[0] / 12.0) * row[4] for row in cashflows) / balance
