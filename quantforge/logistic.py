"""Logistic regression by iteratively reweighted least squares (IRLS).

Models a binary outcome ``y in {0, 1}`` as ``P(y=1|x) = sigmoid(x . beta)``. The
maximum-likelihood coefficients are found by Newton-Raphson, which for the
logistic likelihood is IRLS: each step solves a weighted least-squares problem
with weights ``p(1-p)``. A ridge term stabilizes the Hessian against perfect
separation (where the unpenalized MLE diverges). Used for default-probability and
other classification models. Pure standard library.
"""

import math

from .portopt import _invert


def _sigmoid(z):
    if z >= 0.0:
        return 1.0 / (1.0 + math.exp(-z))
    ez = math.exp(z)
    return ez / (1.0 + ez)


def fit_logistic(X, y, add_intercept=True, max_iter=100, tol=1e-8, ridge=1e-8):
    """Fit a logistic regression by IRLS / Newton-Raphson.

    Parameters
    ----------
    X : list[list[float]]
        Design matrix, ``n`` rows of features.
    y : list[float]
        Binary outcomes (0 or 1).
    add_intercept : bool
        Prepend an intercept column.
    max_iter, tol : int, float
        Newton iteration cap and convergence tolerance on the coefficient step.
    ridge : float
        Small L2 penalty on the Hessian for numerical stability under separation.

    Returns
    -------
    dict
        ``coefficients`` (intercept first if added), ``iterations``,
        ``log_likelihood``, ``converged``.
    """
    n = len(y)
    if n == 0:
        raise ValueError("need at least one observation")
    if len(X) != n:
        raise ValueError("X and y must have the same number of rows")
    if any(v not in (0, 1, 0.0, 1.0) for v in y):
        raise ValueError("y must be binary (0 or 1)")

    if add_intercept:
        design = [[1.0] + list(map(float, row)) for row in X]
    else:
        design = [list(map(float, row)) for row in X]
    p = len(design[0])
    beta = [0.0] * p

    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        # Predicted probabilities and IRLS weights.
        probs = [_sigmoid(sum(design[t][j] * beta[j] for j in range(p)))
                 for t in range(n)]
        w = [max(pi * (1.0 - pi), 1e-12) for pi in probs]

        # Gradient X'(y - p) and Hessian X'WX (+ ridge).
        grad = [sum(design[t][j] * (y[t] - probs[t]) for t in range(n))
                for j in range(p)]
        hess = [[sum(design[t][i] * w[t] * design[t][j] for t in range(n))
                 for j in range(p)] for i in range(p)]
        for i in range(p):
            hess[i][i] += ridge
        inv = _invert(hess)
        step = [sum(inv[i][j] * grad[j] for j in range(p)) for i in range(p)]
        beta = [beta[i] + step[i] for i in range(p)]
        if max(abs(s) for s in step) < tol:
            converged = True
            break

    # Log-likelihood at the solution.
    ll = 0.0
    for t in range(n):
        z = sum(design[t][j] * beta[j] for j in range(p))
        pt = _sigmoid(z)
        pt = min(max(pt, 1e-15), 1.0 - 1e-15)
        ll += y[t] * math.log(pt) + (1.0 - y[t]) * math.log(1.0 - pt)

    return {
        "coefficients": beta,
        "iterations": it,
        "log_likelihood": ll,
        "converged": converged,
    }


def predict_proba(model, X, add_intercept=True):
    """Predicted P(y=1) for each row of ``X`` under a fitted logistic model."""
    beta = model["coefficients"]
    rows = ([[1.0] + list(map(float, r)) for r in X] if add_intercept
            else [list(map(float, r)) for r in X])
    return [_sigmoid(sum(rows[t][j] * beta[j] for j in range(len(beta))))
            for t in range(len(rows))]
