"""Poisson regression: a generalized linear model for count data.

Linear regression assumes Gaussian, constant-variance errors -- wrong for counts, which
are non-negative integers whose variance grows with the mean. Poisson regression models
``E[y | x] = exp(x' beta)`` (a log link, so the rate is always positive) and fits by
maximum likelihood. The MLE solves by iteratively reweighted least squares: at each step
the working response and weights come from the current fitted rates, and a weighted
least-squares solve updates the coefficients. Standard for event counts, claim
frequencies, and arrival rates. Builds on the library's WLS. Pure standard library.
"""

import math

from .wls import weighted_least_squares


def poisson_regression(X, y, add_intercept=True, max_iter=50, tol=1e-8):
    """Fit a Poisson GLM (log link) by IRLS.

    ``y`` are non-negative counts. Returns a dict with ``coefficients`` (intercept first
    if added), ``n_iter`` and the ``log_likelihood`` at convergence. Predicted rate for
    a row is ``exp(x' beta)``. Uses Fisher-scoring IRLS: working response
    ``z = eta + (y - mu)/mu`` with weights ``mu``.
    """
    n = len(y)
    if n != len(X):
        raise ValueError("X and y must have equal length")
    if n == 0:
        raise ValueError("need at least one observation")
    if any(v < 0 for v in y):
        raise ValueError("counts must be non-negative")

    design = [[1.0] + list(map(float, row)) for row in X] if add_intercept \
        else [list(map(float, row)) for row in X]
    p = len(design[0])
    # Start intercept at log(mean count), other slopes at 0.
    ybar = sum(y) / n
    beta = [0.0] * p
    beta[0] = math.log(ybar) if (add_intercept and ybar > 0) else 0.0

    for it in range(1, max_iter + 1):
        eta = [sum(design[t][j] * beta[j] for j in range(p)) for t in range(n)]
        mu = [math.exp(min(e, 700.0)) for e in eta]         # guard overflow
        # Working response and IRLS weights.
        z = [eta[t] + (y[t] - mu[t]) / mu[t] if mu[t] > 0 else eta[t] for t in range(n)]
        w = [mu[t] for t in range(n)]
        # Weighted LS of z on the design (design already includes intercept).
        fit = weighted_least_squares(design, z, w, add_intercept=False)
        new_beta = fit["coefficients"]
        change = max(abs(new_beta[j] - beta[j]) for j in range(p))
        beta = new_beta
        if change < tol:
            break

    # Log-likelihood (drop the constant log(y!) term; sufficient for comparison).
    eta = [sum(design[t][j] * beta[j] for j in range(p)) for t in range(n)]
    ll = sum(y[t] * eta[t] - math.exp(min(eta[t], 700.0)) for t in range(n))
    return {"coefficients": beta, "n_iter": it, "log_likelihood": ll}


def poisson_predict(model, X_query, add_intercept=True):
    """Predicted rates ``exp(x' beta)`` for rows ``X_query``."""
    beta = model["coefficients"]
    rows = [[1.0] + list(map(float, r)) for r in X_query] if add_intercept \
        else [list(map(float, r)) for r in X_query]
    p = len(beta)
    return [math.exp(min(sum(row[j] * beta[j] for j in range(p)), 700.0)) for row in rows]
