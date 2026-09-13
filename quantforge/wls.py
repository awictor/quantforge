"""Weighted and generalized least squares.

OLS treats every observation equally; when observations differ in reliability or the
errors are correlated, weighting them correctly is the efficient (minimum-variance)
estimator:

  * ``weighted_least_squares`` -- minimizes ``sum w_t (y_t - x_t beta)^2``. Weights
    proportional to ``1 / var(eps_t)`` down-weight noisy observations; equal weights
    reproduce OLS exactly.
  * ``generalized_least_squares`` -- for a known error covariance ``Sigma``, solves
    ``beta = (X' Sigma^{-1} X)^{-1} X' Sigma^{-1} y`` by whitening the system with the
    Cholesky factor of ``Sigma``, then running OLS on the transformed data.

Both return coefficients and their standard errors. Pure standard library; builds on
the library's matrix inverse and Cholesky.
"""

import math

from .portopt import _invert
from .linalg import cholesky


def _design(X, add_intercept):
    if add_intercept:
        return [[1.0] + list(map(float, row)) for row in X]
    return [list(map(float, row)) for row in X]


def weighted_least_squares(X, y, weights, add_intercept=True):
    """Weighted least squares: minimize ``sum w_t (y_t - x_t beta)^2``.

    ``weights`` is a length-``n`` list of non-negative weights (larger = more trusted).
    Returns a dict with ``coefficients``, ``std_errors`` (using the weighted residual
    variance), ``residuals`` and ``r_squared`` (weighted). Equal weights reproduce OLS.
    """
    n = len(y)
    if n != len(X) or n != len(weights):
        raise ValueError("X, y, weights must have the same number of rows")
    if any(w < 0 for w in weights):
        raise ValueError("weights must be non-negative")
    design = _design(X, add_intercept)
    p = len(design[0])
    if n <= p:
        raise ValueError("need more observations than parameters")

    # Weighted normal equations: (X' W X) beta = X' W y.
    xtwx = [[sum(weights[t] * design[t][i] * design[t][j] for t in range(n))
             for j in range(p)] for i in range(p)]
    xtwy = [sum(weights[t] * design[t][i] * y[t] for t in range(n)) for i in range(p)]
    inv = _invert(xtwx)
    beta = [sum(inv[i][j] * xtwy[j] for j in range(p)) for i in range(p)]

    resid = [y[t] - sum(design[t][j] * beta[j] for j in range(p)) for t in range(n)]
    dof = n - p
    # Weighted residual variance.
    s2 = sum(weights[t] * resid[t] ** 2 for t in range(n)) / dof
    se = [math.sqrt(s2 * inv[i][i]) if inv[i][i] > 0 else 0.0 for i in range(p)]

    wmean = sum(weights[t] * y[t] for t in range(n)) / sum(weights)
    ss_tot = sum(weights[t] * (y[t] - wmean) ** 2 for t in range(n))
    ss_res = sum(weights[t] * resid[t] ** 2 for t in range(n))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return {"coefficients": beta, "std_errors": se, "residuals": resid,
            "r_squared": r2}


def _forward_solve(L, b):
    """Solve L z = b for lower-triangular L."""
    n = len(b)
    z = [0.0] * n
    for i in range(n):
        z[i] = (b[i] - sum(L[i][j] * z[j] for j in range(i))) / L[i][i]
    return z


def generalized_least_squares(X, y, cov, add_intercept=True):
    """Generalized least squares for a known error covariance ``cov`` (Sigma).

    Whitens the system with the Cholesky factor of ``Sigma`` (``Sigma = L L'``), so
    ``L^{-1} y = L^{-1} X beta + white noise``, then applies OLS. Returns a dict with
    ``coefficients`` and ``std_errors``. A diagonal ``cov`` reproduces weighted least
    squares with ``weights = 1 / diag(cov)``.
    """
    n = len(y)
    if n != len(X) or n != len(cov):
        raise ValueError("X, y, cov must have matching dimensions")
    design = _design(X, add_intercept)
    p = len(design[0])
    if n <= p:
        raise ValueError("need more observations than parameters")

    L = cholesky(cov)                        # Sigma = L L'
    # Whiten: y* = L^{-1} y, X* = L^{-1} X (column by column).
    y_w = _forward_solve(L, y)
    x_w = []
    for j in range(p):
        col = [design[t][j] for t in range(n)]
        x_w.append(_forward_solve(L, col))
    # Rebuild whitened design as rows.
    xw = [[x_w[j][t] for j in range(p)] for t in range(n)]

    xtx = [[sum(xw[t][i] * xw[t][j] for t in range(n)) for j in range(p)]
           for i in range(p)]
    xty = [sum(xw[t][i] * y_w[t] for t in range(n)) for i in range(p)]
    inv = _invert(xtx)
    beta = [sum(inv[i][j] * xty[j] for j in range(p)) for i in range(p)]

    resid_w = [y_w[t] - sum(xw[t][j] * beta[j] for j in range(p)) for t in range(n)]
    s2 = sum(r * r for r in resid_w) / (n - p)
    se = [math.sqrt(s2 * inv[i][i]) if inv[i][i] > 0 else 0.0 for i in range(p)]
    return {"coefficients": beta, "std_errors": se}
