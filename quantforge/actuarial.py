"""Life-contingent actuarial functions: survival, annuities, and insurance.

Given one-year survival probabilities ``p_x`` from a life table, the ``k``-year
survival is ``kp_x = prod_{j<k} p_{x+j}``. Standard expected present values at a
flat annual interest rate ``i`` (discount ``v = 1/(1+i)``) are

    life annuity-due   a-due_x = sum_k v^k * kp_x
    term insurance     A_term  = sum_k v^{k+1} * kp_x * q_{x+k}
    whole life         A_x     = term over the whole table
    endowment          A_endow = term insurance + pure endowment

These satisfy the identity ``A_x = 1 - d * a-due_x`` with ``d = i v`` for a whole-
life contract over the full table. Pure standard library.
"""

import math


def survival_probabilities(one_year_survival):
    """Cumulative survival ``[0p_x, 1p_x, 2p_x, ...]`` from one-year ``p_x`` values.

    ``kp_x = prod_{j<k} p_{x+j}``, starting at ``0p_x = 1``. The returned list has
    one more entry than the input (the leading 1). Survival is non-increasing.
    """
    if any(not (0.0 <= p <= 1.0) for p in one_year_survival):
        raise ValueError("survival probabilities must be in [0, 1]")
    cum = [1.0]
    acc = 1.0
    for p in one_year_survival:
        acc *= p
        cum.append(acc)
    return cum


def life_annuity_due(one_year_survival, i):
    """Expected present value of a unit life annuity-due.

    ``a-due = sum_k v^k * kp_x`` paying 1 at the start of each year while alive,
    ``v = 1/(1+i)``. Falls as interest or mortality rises.
    """
    if i <= -1.0:
        raise ValueError("interest rate must exceed -100%")
    v = 1.0 / (1.0 + i)
    cum = survival_probabilities(one_year_survival)
    return sum(v ** k * cum[k] for k in range(len(cum)))


def term_insurance(one_year_survival, i, term=None):
    """EPV of a unit term insurance paying 1 at the end of the year of death.

    ``A = sum_k v^{k+1} * kp_x * q_{x+k}`` over the first ``term`` years (default:
    the whole table). ``q = 1 - p`` is the one-year death probability.
    """
    if i <= -1.0:
        raise ValueError("interest rate must exceed -100%")
    v = 1.0 / (1.0 + i)
    cum = survival_probabilities(one_year_survival)
    n = len(one_year_survival) if term is None else min(term, len(one_year_survival))
    total = 0.0
    for k in range(n):
        q = 1.0 - one_year_survival[k]
        total += v ** (k + 1) * cum[k] * q
    return total


def whole_life_insurance(one_year_survival, i):
    """EPV of whole-life insurance: :func:`term_insurance` over the whole table."""
    return term_insurance(one_year_survival, i)


def pure_endowment(one_year_survival, i, term):
    """EPV of a unit pure endowment: ``v^n * np_x`` (pays 1 iff alive at ``n``)."""
    if i <= -1.0:
        raise ValueError("interest rate must exceed -100%")
    v = 1.0 / (1.0 + i)
    cum = survival_probabilities(one_year_survival)
    n = min(term, len(cum) - 1)
    return v ** n * cum[n]


def endowment_insurance(one_year_survival, i, term):
    """EPV of an endowment: term insurance plus a pure endowment at ``term``."""
    return (term_insurance(one_year_survival, i, term)
            + pure_endowment(one_year_survival, i, term))
