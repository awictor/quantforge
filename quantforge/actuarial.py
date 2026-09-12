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


def cat_layer_loss(gross_loss, attachment, exhaustion):
    """Loss ceded to a reinsurance / cat-bond layer ``[attachment, exhaustion]``.

    ``min(max(gross_loss - attachment, 0), exhaustion - attachment)`` -- zero below
    the attachment point, rising one-for-one through the layer, capped at the layer
    width above exhaustion.
    """
    if attachment < 0 or exhaustion <= attachment:
        raise ValueError("require 0 <= attachment < exhaustion")
    return min(max(gross_loss - attachment, 0.0), exhaustion - attachment)


def cat_expected_loss(loss_scenarios, attachment, exhaustion):
    """Expected layer loss over equally-likely loss scenarios (as a fraction).

    Averages :func:`cat_layer_loss` across ``loss_scenarios`` and divides by the
    layer width, giving the expected loss as a fraction of the layer notional in
    ``[0, 1]`` -- the cat bond's expected loss rate.
    """
    if not loss_scenarios:
        raise ValueError("need at least one loss scenario")
    width = exhaustion - attachment
    avg = sum(cat_layer_loss(x, attachment, exhaustion)
              for x in loss_scenarios) / len(loss_scenarios)
    return avg / width


def cat_bond_spread(expected_loss_rate, risk_load=1.0):
    """Fair coupon spread of a cat bond: expected loss rate times a risk load.

    Investors demand a spread above the expected loss to bear the (undiversifiable,
    fat-tailed) catastrophe risk: ``spread = (1 + risk_load) * expected_loss_rate``
    (``risk_load = 0`` is the actuarially fair spread). At or above the expected
    loss rate.
    """
    if not (0.0 <= expected_loss_rate <= 1.0):
        raise ValueError("expected_loss_rate must be in [0, 1]")
    if risk_load < 0:
        raise ValueError("risk_load must be non-negative")
    return (1.0 + risk_load) * expected_loss_rate


def cat_bond_price(principal, coupon_rate, expected_loss_rate, r, maturity,
                   risk_free_spread=0.0):
    """Present value of a single-period cat bond.

    Pays ``coupon_rate`` on the principal and returns the principal at maturity
    unless a triggering event erodes it by the expected loss. Discounts the
    expected principal repayment ``principal * (1 - expected_loss_rate)`` and the
    coupon at ``r + risk_free_spread``. Falls as the expected loss rises.
    """
    if not (0.0 <= expected_loss_rate <= 1.0):
        raise ValueError("expected_loss_rate must be in [0, 1]")
    if maturity < 0:
        raise ValueError("maturity must be non-negative")
    disc = math.exp(-(r + risk_free_spread) * maturity)
    coupon = principal * coupon_rate * maturity
    repay = principal * (1.0 - expected_loss_rate)
    return disc * (coupon + repay)


def gompertz_makeham_hazard(age, a, b, c):
    """Gompertz-Makeham force of mortality ``mu(x) = a + b * c^x``.

    ``a`` is the age-independent (accident) component and ``b c^x`` the
    exponentially-rising Gompertz term. Increasing in age for ``c > 1``.
    """
    if a < 0 or b < 0:
        raise ValueError("a and b must be non-negative")
    if c <= 0:
        raise ValueError("c must be positive")
    return a + b * c ** age


def gompertz_makeham_survival(age, years, a, b, c):
    """Survival probability over ``years`` under Gompertz-Makeham mortality.

    Integrates the force of mortality from ``age`` to ``age + years``:

        tp_x = exp(-a t - (b / ln c) c^x (c^t - 1)),   t = years

    (the closed-form integral of ``a + b c^s``). Falls monotonically with the
    horizon; the ``c -> 1`` limit uses the exponential (Makeham-only) form.
    """
    if years < 0:
        raise ValueError("years must be non-negative")
    if a < 0 or b < 0:
        raise ValueError("a and b must be non-negative")
    if c <= 0:
        raise ValueError("c must be positive")
    makeham = a * years
    if abs(c - 1.0) < 1e-12:
        gompertz = b * years
    else:
        gompertz = b / math.log(c) * c ** age * (c ** years - 1.0)
    return math.exp(-(makeham + gompertz))


def gompertz_makeham_survival_curve(age, n_years, a, b, c):
    """One-year survival probabilities ``[p_x, p_{x+1}, ...]`` for ``n_years``.

    Each entry is the one-year Gompertz-Makeham survival at successive ages, ready
    to feed the life-table functions (:func:`life_annuity_due`, etc.).
    """
    if n_years < 1:
        raise ValueError("n_years must be a positive integer")
    return [gompertz_makeham_survival(age + k, 1.0, a, b, c) for k in range(n_years)]


def curtate_life_expectancy(one_year_survival):
    """Curtate expectation of life ``e_x = sum_{k>=1} kp_x`` (whole years).

    The expected number of complete future years lived, the sum of the cumulative
    survival probabilities beyond time zero.
    """
    cum = survival_probabilities(one_year_survival)
    return sum(cum[1:])


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


def temporary_life_annuity_due(one_year_survival, i, term):
    """EPV of an ``n``-year temporary life annuity-due.

    ``a-due_{x:n} = sum_{k<n} v^k * kp_x`` -- pays 1 at the start of each year while
    alive, for at most ``term`` years. Below the whole-life
    :func:`life_annuity_due` and rising to it as ``term`` grows.
    """
    if i <= -1.0:
        raise ValueError("interest rate must exceed -100%")
    v = 1.0 / (1.0 + i)
    cum = survival_probabilities(one_year_survival)
    n = min(term, len(cum))
    return sum(v ** k * cum[k] for k in range(n))


def net_level_premium(one_year_survival, i, term=None):
    """Net annual premium for a (term or whole-life) unit insurance.

    By the equivalence principle the level premium equates the EPV of premiums
    (a life annuity-due) to the EPV of benefits (the insurance):

        P = A / a-due

    Uses :func:`whole_life_insurance` over ``a-due`` for whole life (``term`` None)
    or :func:`endowment_insurance` over the temporary annuity for an ``term``-year
    endowment. The premium the insurer must charge to break even.
    """
    if term is None:
        benefit = whole_life_insurance(one_year_survival, i)
        annuity = life_annuity_due(one_year_survival, i)
    else:
        benefit = endowment_insurance(one_year_survival, i, term)
        annuity = temporary_life_annuity_due(one_year_survival, i, term)
    if annuity <= 0.0:
        raise ValueError("annuity EPV must be positive")
    return benefit / annuity
