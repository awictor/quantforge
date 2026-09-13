"""Robust OLS standard errors: White (heteroskedasticity) and Newey-West (HAC).

The textbook OLS standard errors assume the errors are homoskedastic and serially
uncorrelated. When they are not -- volatility clustering, autocorrelated residuals --
those errors are wrong and t-statistics mislead. The sandwich estimator fixes it:

    Var(beta) = (X'X)^{-1} S (X'X)^{-1},

with the meat ``S`` capturing the true error structure. ``white_hc0`` uses
``S = sum e_t^2 x_t x_t'`` (heteroskedasticity-consistent); ``newey_west`` adds
Bartlett-weighted autocovariance terms out to ``lags`` for autocorrelation. Returns the
robust standard errors and t-statistics for each coefficient. Pure standard library;
builds on the OLS normal-equations solve.
"""

import math

from .portopt import _invert


def _design(X, add_intercept):
    if add_intercept:
        return [[1.0] + list(map(float, row)) for row in X]
    return [list(map(float, row)) for row in X]


def _fit(design, y):
    n = len(y)
    p = len(design[0])
    xtx = [[sum(design[t][i] * design[t][j] for t in range(n)) for j in range(p)]
           for i in range(p)]
    xty = [sum(design[t][i] * y[t] for t in range(n)) for i in range(p)]
    inv = _invert(xtx)
    beta = [sum(inv[i][j] * xty[j] for j in range(p)) for i in range(p)]
    resid = [y[t] - sum(design[t][j] * beta[j] for j in range(p)) for t in range(n)]
    return beta, resid, inv, n, p


def _sandwich(design, inv, meat, n, p):
    # (X'X)^{-1} S (X'X)^{-1}.
    tmp = [[sum(inv[i][k] * meat[k][j] for k in range(p)) for j in range(p)]
           for i in range(p)]
    cov = [[sum(tmp[i][k] * inv[k][j] for k in range(p)) for j in range(p)]
           for i in range(p)]
    return cov


def white_hc0(X, y, add_intercept=True):
    """White (HC0) heteroskedasticity-consistent OLS standard errors.

    Returns a dict with ``coefficients``, robust ``std_errors``, ``t_stats`` and the
    full ``cov`` matrix. Valid when errors are heteroskedastic but not autocorrelated.
    """
    design = _design(X, add_intercept)
    beta, resid, inv, n, p = _fit(design, y)
    meat = [[sum(resid[t] ** 2 * design[t][i] * design[t][j] for t in range(n))
             for j in range(p)] for i in range(p)]
    cov = _sandwich(design, inv, meat, n, p)
    se = [math.sqrt(cov[i][i]) if cov[i][i] > 0 else 0.0 for i in range(p)]
    tstat = [beta[i] / se[i] if se[i] > 0 else 0.0 for i in range(p)]
    return {"coefficients": beta, "std_errors": se, "t_stats": tstat, "cov": cov}


def newey_west(X, y, lags, add_intercept=True):
    """Newey-West HAC OLS standard errors (heteroskedasticity + autocorrelation).

    Adds Bartlett-weighted cross-products of the score vectors ``e_t x_t`` out to
    ``lags`` lags to the White meat, so the covariance is consistent under both
    heteroskedasticity and serial correlation. Returns the same dict shape as
    :func:`white_hc0`. ``lags = 0`` reduces exactly to White (HC0).
    """
    if lags < 0:
        raise ValueError("lags must be non-negative")
    design = _design(X, add_intercept)
    beta, resid, inv, n, p = _fit(design, y)
    if lags >= n:
        raise ValueError("lags must be less than the sample size")

    # Score vectors g_t = e_t * x_t.
    g = [[resid[t] * design[t][i] for i in range(p)] for t in range(n)]

    # Gamma_0.
    meat = [[sum(g[t][i] * g[t][j] for t in range(n)) for j in range(p)]
            for i in range(p)]
    # Bartlett-weighted lags.
    for L in range(1, lags + 1):
        w = 1.0 - L / (lags + 1.0)
        for i in range(p):
            for j in range(p):
                gamma = sum(g[t][i] * g[t - L][j] for t in range(L, n))
                meat[i][j] += w * (gamma + sum(g[t - L][i] * g[t][j]
                                               for t in range(L, n)))
    cov = _sandwich(design, inv, meat, n, p)
    se = [math.sqrt(cov[i][i]) if cov[i][i] > 0 else 0.0 for i in range(p)]
    tstat = [beta[i] / se[i] if se[i] > 0 else 0.0 for i in range(p)]
    return {"coefficients": beta, "std_errors": se, "t_stats": tstat, "cov": cov}
