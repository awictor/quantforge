"""Ridge (Tikhonov) regression: L2-penalized least squares.

Ridge adds an L2 penalty to the OLS objective,

    minimize  ||y - X beta||^2 + lambda ||beta||^2,

giving the closed form ``beta = (X'X + lambda I)^{-1} X'y``. The penalty shrinks
the coefficients toward zero, trading a little bias for much lower variance -- and
it makes the system solvable even when ``X'X`` is singular (perfectly collinear
regressors), where plain OLS fails. The intercept is added by default and is *not*
penalized (only the slopes are shrunk). ``lambda = 0`` reproduces OLS. Pure
standard library.
"""

from .portopt import _invert


def ridge_regression(X, y, alpha=1.0, add_intercept=True):
    """Fit an L2-penalized (ridge) regression.

    Parameters
    ----------
    X : list[list[float]]
        Design matrix, ``n`` rows of regressors.
    y : list[float]
        Response vector.
    alpha : float
        Ridge penalty ``lambda`` (>= 0). 0 reproduces OLS; larger shrinks the
        slope coefficients toward zero.
    add_intercept : bool
        Prepend an (unpenalized) intercept column.

    Returns
    -------
    dict
        ``coefficients`` (intercept first if added), ``fitted``, ``residuals``,
        ``r_squared``.
    """
    n = len(y)
    if n == 0:
        raise ValueError("need at least one observation")
    if len(X) != n:
        raise ValueError("X and y must have the same number of rows")
    if alpha < 0:
        raise ValueError("alpha must be non-negative")

    if add_intercept:
        design = [[1.0] + list(map(float, row)) for row in X]
    else:
        design = [list(map(float, row)) for row in X]
    p = len(design[0])

    xtx = [[sum(design[t][i] * design[t][j] for t in range(n)) for j in range(p)]
           for i in range(p)]
    # Add lambda to the diagonal, skipping the intercept term if present.
    start = 1 if add_intercept else 0
    for i in range(start, p):
        xtx[i][i] += alpha
    xty = [sum(design[t][i] * y[t] for t in range(n)) for i in range(p)]
    inv = _invert(xtx)
    beta = [sum(inv[i][j] * xty[j] for j in range(p)) for i in range(p)]

    fitted = [sum(design[t][i] * beta[i] for i in range(p)) for t in range(n)]
    resid = [y[t] - fitted[t] for t in range(n)]
    ss_res = sum(r * r for r in resid)
    ybar = sum(y) / n
    ss_tot = sum((yt - ybar) ** 2 for yt in y)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return {
        "coefficients": beta,
        "fitted": fitted,
        "residuals": resid,
        "r_squared": r2,
    }
