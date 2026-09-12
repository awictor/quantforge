"""Ordinary least squares with regression diagnostics.

Fits ``y = X beta + eps`` by the normal equations ``beta = (X'X)^{-1} X'y`` and
reports the standard diagnostics: coefficient standard errors, t-statistics,
R-squared and adjusted R-squared, and the overall F-statistic. An intercept column
is added by default. Pure standard library (uses a Gauss-Jordan inverse).
"""

import math

from .portopt import _invert


def ols_fit(X, y, add_intercept=True):
    """Fit an OLS regression and return coefficients with diagnostics.

    Parameters
    ----------
    X : list[list[float]]
        Design matrix, ``n`` rows of ``k`` regressors (no intercept column unless
        ``add_intercept=False`` and you supply your own).
    y : list[float]
        Response vector of length ``n``.
    add_intercept : bool
        Prepend a column of ones (the default).

    Returns
    -------
    dict
        ``coefficients`` (intercept first if added), ``std_errors``, ``t_stats``,
        ``r_squared``, ``adj_r_squared``, ``f_stat``, ``residuals``, ``n_obs``,
        ``df_resid``.
    """
    n = len(y)
    if n == 0:
        raise ValueError("need at least one observation")
    if len(X) != n:
        raise ValueError("X and y must have the same number of rows")

    if add_intercept:
        design = [[1.0] + list(map(float, row)) for row in X]
    else:
        design = [list(map(float, row)) for row in X]
    p = len(design[0])
    if n <= p:
        raise ValueError("need more observations than parameters")

    # Normal equations.
    xtx = [[sum(design[t][i] * design[t][j] for t in range(n)) for j in range(p)]
           for i in range(p)]
    xty = [sum(design[t][i] * y[t] for t in range(n)) for i in range(p)]
    inv = _invert(xtx)
    beta = [sum(inv[i][j] * xty[j] for j in range(p)) for i in range(p)]

    fitted = [sum(design[t][i] * beta[i] for i in range(p)) for t in range(n)]
    resid = [y[t] - fitted[t] for t in range(n)]
    ss_res = sum(r * r for r in resid)
    ybar = sum(y) / n
    ss_tot = sum((yt - ybar) ** 2 for yt in y)

    df_resid = n - p
    sigma2 = ss_res / df_resid if df_resid > 0 else 0.0
    std_err = [math.sqrt(sigma2 * inv[i][i]) if inv[i][i] > 0 else 0.0
               for i in range(p)]
    t_stats = [beta[i] / std_err[i] if std_err[i] > 0 else 0.0 for i in range(p)]

    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    # Adjusted R^2 and F use the number of non-intercept regressors.
    k = p - 1 if add_intercept else p
    adj_r2 = (1.0 - (1.0 - r2) * (n - 1) / df_resid) if df_resid > 0 else r2
    if k > 0 and ss_res > 0:
        f_stat = ((ss_tot - ss_res) / k) / (ss_res / df_resid)
    else:
        f_stat = float("inf") if ss_res == 0 else 0.0

    return {
        "coefficients": beta,
        "std_errors": std_err,
        "t_stats": t_stats,
        "r_squared": r2,
        "adj_r_squared": adj_r2,
        "f_stat": f_stat,
        "residuals": resid,
        "n_obs": n,
        "df_resid": df_resid,
    }
