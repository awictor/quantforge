"""Multi-factor (OLS) return regression: alpha, factor betas, and fit.

Regresses an asset's excess returns on one or more factor return series
(market, size, value, momentum, ...) by ordinary least squares, returning the
intercept (alpha), the factor loadings (betas), the R-squared, and the residual
volatility. Solved via the normal equations with a pure-Python Gaussian
elimination -- no external linear-algebra dependency.
"""

import math


def _solve(matrix, rhs):
    """Solve a small dense linear system by Gaussian elimination with pivoting."""
    n = len(matrix)
    a = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[piv][col]) < 1e-15:
            raise ValueError("singular system (collinear factors?)")
        a[col], a[piv] = a[piv], a[col]
        pivot = a[col][col]
        for j in range(col, n + 1):
            a[col][j] /= pivot
        for r in range(n):
            if r != col:
                factor = a[r][col]
                for j in range(col, n + 1):
                    a[r][j] -= factor * a[col][j]
    return [a[i][n] for i in range(n)]


def factor_regression(asset_returns, factor_returns):
    """OLS regression of ``asset_returns`` on ``factor_returns`` (list of series).

    ``factor_returns`` is a list of equal-length factor return series. Returns a
    dict with ``alpha`` (intercept), ``betas`` (one loading per factor),
    ``r_squared``, and ``residual_vol`` (standard deviation of the residuals).
    Fits ``r_t = alpha + sum_k beta_k f_{k,t} + eps_t`` by minimizing the squared
    residuals.
    """
    n = len(asset_returns)
    if n == 0:
        raise ValueError("need at least one observation")
    k = len(factor_returns)
    for f in factor_returns:
        if len(f) != n:
            raise ValueError("all series must have equal length")
    # Design matrix columns: intercept + each factor.
    cols = [[1.0] * n] + [list(f) for f in factor_returns]
    p = k + 1
    if n < p:
        raise ValueError("need at least as many observations as parameters")
    # Normal equations X'X beta = X'y.
    xtx = [[sum(cols[i][t] * cols[j][t] for t in range(n)) for j in range(p)]
           for i in range(p)]
    xty = [sum(cols[i][t] * asset_returns[t] for t in range(n)) for i in range(p)]
    coeffs = _solve(xtx, xty)
    alpha = coeffs[0]
    betas = coeffs[1:]
    # Fit statistics.
    mean_y = sum(asset_returns) / n
    ss_tot = sum((y - mean_y) ** 2 for y in asset_returns)
    resid = []
    for t in range(n):
        pred = alpha + sum(betas[i] * factor_returns[i][t] for i in range(k))
        resid.append(asset_returns[t] - pred)
    ss_res = sum(e * e for e in resid)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    residual_vol = math.sqrt(ss_res / n)
    return {"alpha": alpha, "betas": betas, "r_squared": r_squared,
            "residual_vol": residual_vol}


def factor_expected_return(alpha, betas, factor_premia):
    """Expected return from a fitted factor model: ``alpha + sum beta_k premium_k``."""
    if len(betas) != len(factor_premia):
        raise ValueError("betas and factor_premia must align")
    return alpha + sum(b * p for b, p in zip(betas, factor_premia))
