"""Bayesian linear regression with a conjugate Gaussian prior.

Instead of a single least-squares point estimate, Bayesian regression returns a
*distribution* over the coefficients. With a Gaussian prior ``beta ~ N(0, alpha^{-1} I)``
and known noise precision ``beta_noise``, the posterior is Gaussian in closed form:

    S = (alpha I + beta_noise X'X)^{-1},   mean = beta_noise S X'y.

The posterior mean coincides with ridge regression (``lambda = alpha / beta_noise``),
but the posterior *covariance* also gives coefficient credible intervals and, for a new
point, a predictive variance -- the model's own uncertainty, which a point estimate
can't express. Builds on the library's matrix inverse. Pure standard library.
"""

import math

from .portopt import _invert


def bayesian_linear_regression(X, y, alpha=1.0, beta_noise=1.0, add_intercept=True):
    """Conjugate Bayesian linear regression posterior.

    ``alpha`` is the prior precision (larger = stronger shrinkage toward zero),
    ``beta_noise`` the noise precision (``1 / variance``). Returns a dict with the
    posterior ``mean`` (coefficient vector, intercept first if added), the posterior
    ``covariance`` matrix, and ``std`` (per-coefficient posterior standard deviations).
    The mean equals ridge regression with ``lambda = alpha / beta_noise`` that also
    penalizes the intercept (the prior shrinks *every* coefficient toward zero).
    """
    n = len(y)
    if n != len(X):
        raise ValueError("X and y must have equal length")
    if n == 0:
        raise ValueError("need at least one observation")
    if alpha <= 0 or beta_noise <= 0:
        raise ValueError("alpha and beta_noise must be positive")

    design = [[1.0] + list(map(float, row)) for row in X] if add_intercept \
        else [list(map(float, row)) for row in X]
    p = len(design[0])

    # Posterior precision A = alpha I + beta_noise X'X.
    xtx = [[sum(design[t][i] * design[t][j] for t in range(n)) for j in range(p)]
           for i in range(p)]
    A = [[(alpha if i == j else 0.0) + beta_noise * xtx[i][j] for j in range(p)]
         for i in range(p)]
    S = _invert(A)                                   # posterior covariance
    xty = [sum(design[t][i] * y[t] for t in range(n)) for i in range(p)]
    mean = [beta_noise * sum(S[i][j] * xty[j] for j in range(p)) for i in range(p)]
    std = [math.sqrt(S[i][i]) if S[i][i] > 0 else 0.0 for i in range(p)]
    return {"mean": mean, "covariance": S, "std": std,
            "_design_p": p, "_add_intercept": add_intercept, "_beta_noise": beta_noise}


def bayesian_predict(model, x_row):
    """Predictive mean and variance for a new point ``x_row``.

    Returns ``(mean, variance)`` where the variance is the *predictive* variance
    ``1/beta_noise + x' S x`` -- observation noise plus the posterior uncertainty in
    the coefficients (so it widens where data is sparse). ``x_row`` excludes the
    intercept if the model was fit with one.
    """
    S = model["covariance"]
    mean_beta = model["mean"]
    p = model["_design_p"]
    xv = ([1.0] + list(map(float, x_row))) if model["_add_intercept"] \
        else list(map(float, x_row))
    if len(xv) != p:
        raise ValueError("x_row has wrong dimension")
    mean = sum(mean_beta[i] * xv[i] for i in range(p))
    var = 1.0 / model["_beta_noise"] + sum(
        xv[i] * S[i][j] * xv[j] for i in range(p) for j in range(p))
    return mean, var
