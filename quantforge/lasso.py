"""LASSO regression by cyclic coordinate descent.

The LASSO adds an L1 penalty to least squares,

    min_beta (1/2n) ||y - X beta||^2 + alpha ||beta||_1,

which shrinks coefficients and drives some exactly to zero -- so it both regularizes
and selects features, unlike ridge's L2 penalty. This module solves it by cyclic
coordinate descent with soft-thresholding on standardized features; the intercept is
not penalized. Coefficients are returned on the original (unstandardized) scale.
Pure standard library.
"""


def _soft_threshold(z, gamma):
    if z > gamma:
        return z - gamma
    if z < -gamma:
        return z + gamma
    return 0.0


def lasso_regression(X, y, alpha=1.0, max_iter=1000, tol=1e-8):
    """Fit a LASSO regression ``y ~ X beta`` by coordinate descent.

    Returns ``[intercept, b_1, ..., b_p]`` on the original feature scale. ``alpha``
    is the L1 penalty strength: ``alpha = 0`` recovers ordinary least squares, and a
    large ``alpha`` drives all slopes to zero (the fit collapses to the mean of
    ``y``). Features are standardized internally so the penalty applies evenly; the
    intercept is never penalized. A slope set to exactly zero has been selected out.
    """
    n = len(y)
    if n == 0 or len(X) != n:
        raise ValueError("X and y must have the same number of rows")
    p = len(X[0])
    if any(len(row) != p for row in X):
        raise ValueError("all rows of X must have the same length")
    if alpha < 0.0:
        raise ValueError("alpha must be non-negative")

    # Standardize features (mean 0, unit variance); center y.
    col_mean = [sum(X[i][j] for i in range(n)) / n for j in range(p)]
    col_std = []
    for j in range(p):
        var = sum((X[i][j] - col_mean[j]) ** 2 for i in range(n)) / n
        col_std.append(var ** 0.5 if var > 0 else 1.0)
    Z = [[(X[i][j] - col_mean[j]) / col_std[j] for j in range(p)] for i in range(n)]
    y_mean = sum(y) / n
    yc = [y[i] - y_mean for i in range(n)]

    beta = [0.0] * p
    # Precompute column norms (standardized => each is n).
    for _ in range(max_iter):
        max_change = 0.0
        for j in range(p):
            # Partial residual excluding feature j.
            rho = 0.0
            for i in range(n):
                pred_wo_j = sum(Z[i][k] * beta[k] for k in range(p)) - Z[i][j] * beta[j]
                rho += Z[i][j] * (yc[i] - pred_wo_j)
            new_bj = _soft_threshold(rho / n, alpha) / (1.0)  # col norm/n = 1
            max_change = max(max_change, abs(new_bj - beta[j]))
            beta[j] = new_bj
        if max_change < tol:
            break

    # Map back to the original scale: b_orig_j = beta_j / std_j; adjust intercept.
    slopes = [beta[j] / col_std[j] for j in range(p)]
    intercept = y_mean - sum(slopes[j] * col_mean[j] for j in range(p))
    return [intercept] + slopes
