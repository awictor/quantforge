"""Huber robust regression (M-estimator by iteratively reweighted least squares).

OLS squares every residual, so a single outlier dominates the fit; a fully robust
median fit ignores the magnitude of moderate residuals. Huber's loss interpolates: it
is quadratic for small residuals (efficient, like OLS) and linear beyond a threshold
``delta`` (bounded influence, like an absolute loss). Fitting minimizes that loss by
iteratively reweighted least squares -- each step down-weights points whose scaled
residual exceeds ``delta`` by ``delta / |r|``. The result resists outliers while
staying near-OLS-efficient on clean Gaussian data. Builds on the library's WLS. Pure
standard library.
"""

import math

from .wls import weighted_least_squares


def _mad_scale(residuals):
    """Robust residual scale: MAD / 0.6745 (normal-consistent)."""
    med = _median(residuals)
    mad = _median([abs(r - med) for r in residuals])
    return mad / 0.6744897501960817 if mad > 0 else 1.0


def _median(v):
    s = sorted(v)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else 0.5 * (s[mid - 1] + s[mid])


def huber_regression(X, y, delta=1.345, add_intercept=True, max_iter=50, tol=1e-8):
    """Huber robust regression by IRLS.

    ``delta`` is the residual threshold (in robust-scale units) between the quadratic
    and linear regions; the default 1.345 gives ~95% efficiency at the normal. Returns
    a dict with ``coefficients``, ``n_iter`` and ``scale`` (the final robust residual
    scale). Reduces toward OLS as ``delta`` grows.
    """
    n = len(y)
    if n != len(X):
        raise ValueError("X and y must have equal length")
    if n == 0:
        raise ValueError("need at least one observation")
    if delta <= 0:
        raise ValueError("delta must be positive")

    design = [[1.0] + list(map(float, row)) for row in X] if add_intercept \
        else [list(map(float, row)) for row in X]
    p = len(design[0])
    beta = [0.0] * p
    scale = 1.0

    for it in range(1, max_iter + 1):
        resid = [y[t] - sum(design[t][j] * beta[j] for j in range(p)) for t in range(n)]
        scale = _mad_scale(resid)
        if scale <= 0:
            break
        # Huber weights: 1 for |r/scale| <= delta, else delta*scale/|r|.
        weights = []
        for r in resid:
            a = abs(r) / scale
            weights.append(1.0 if a <= delta else delta / a)
        fit = weighted_least_squares(X, y, weights, add_intercept=add_intercept)
        new_beta = fit["coefficients"]
        change = max(abs(new_beta[j] - beta[j]) for j in range(p))
        beta = new_beta
        if change < tol:
            return {"coefficients": beta, "n_iter": it, "scale": scale}
    return {"coefficients": beta, "n_iter": max_iter, "scale": scale}
