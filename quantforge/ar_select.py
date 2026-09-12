"""AR order selection by information criteria (AIC / BIC).

An AR(p) fit trades goodness-of-fit against parameter count. Fit each candidate
order, evaluate the Gaussian log-likelihood from the estimated innovation
variance, and penalize the number of parameters:

    AIC = -2 logL + 2 k,      BIC = -2 logL + k ln(n),

where ``k = p + 1`` (the AR coefficients plus the noise variance). The order that
minimizes the criterion is selected. BIC's heavier ``ln(n)`` penalty favors more
parsimonious models than AIC. Uses the concentrated Gaussian log-likelihood
``logL = -n/2 (ln(2 pi sigma^2) + 1)``. Pure standard library.
"""

import math

from .ar_model import fit_ar_yule_walker


def _loglik(n, noise_var):
    if noise_var <= 0.0:
        return float("-inf")
    return -0.5 * n * (math.log(2.0 * math.pi * noise_var) + 1.0)


def ar_information_criteria(x, order):
    """AIC and BIC of an AR(``order``) fit to ``x``.

    Returns ``(aic, bic, loglik, noise_variance)`` using the concentrated Gaussian
    log-likelihood and ``k = order + 1`` parameters.
    """
    n = len(x)
    m = fit_ar_yule_walker(x, order)
    ll = _loglik(n, m["noise_variance"])
    k = order + 1
    aic = -2.0 * ll + 2.0 * k
    bic = -2.0 * ll + k * math.log(n)
    return aic, bic, ll, m["noise_variance"]


def select_ar_order(x, max_order=10, criterion="aic"):
    """Select the AR order minimizing AIC or BIC over ``1..max_order``.

    Parameters
    ----------
    x : sequence of float
        The series.
    max_order : int
        Largest candidate order to try.
    criterion : str
        ``"aic"`` or ``"bic"``.

    Returns
    -------
    (best_order, scores) : (int, list[tuple])
        The selected order and a list of ``(order, aic, bic)`` for every candidate.
        BIC tends to pick an order no larger than AIC.
    """
    n = len(x)
    if max_order < 1:
        raise ValueError("max_order must be >= 1")
    if n <= max_order + 1:
        raise ValueError("series too short for max_order")
    if criterion not in ("aic", "bic"):
        raise ValueError("criterion must be 'aic' or 'bic'")

    scores = []
    best_order = 1
    best_val = float("inf")
    for p in range(1, max_order + 1):
        aic, bic, _, _ = ar_information_criteria(x, p)
        scores.append((p, aic, bic))
        val = aic if criterion == "aic" else bic
        if val < best_val:
            best_val = val
            best_order = p
    return best_order, scores
