"""Linear quantile regression by iteratively reweighted least squares.

Ordinary least squares fits the conditional mean; quantile regression fits a
conditional quantile ``tau`` by minimizing the pinball loss

    sum_i rho_tau(y_i - x_i' beta),   rho_tau(u) = u (tau - 1[u < 0]).

This module minimizes it with iteratively-reweighted least squares: the pinball loss
is a weighted absolute loss, and absolute loss is fit by IRLS with weights
``1 / |residual|`` (the asymmetry splits the weight by the sign of the residual). The
fit converges to the quantile-regression solution. An intercept is added by default.
Pure standard library on top of the Gauss-Jordan inverse.
"""

from .portopt import _invert


def quantile_regression(X, y, tau=0.5, add_intercept=True, max_iter=200, tol=1e-8):
    """Fit a linear ``tau``-quantile regression ``y ~ X beta``.

    Returns the coefficient list (intercept first if added). ``tau`` in ``(0, 1)``
    selects the conditional quantile: 0.5 is the median (least-absolute-deviations)
    fit, higher ``tau`` tracks the upper conditional tail. Solved by IRLS on the
    asymmetric absolute loss; a small floor keeps the reweighting stable at zero
    residuals.
    """
    n = len(y)
    if n == 0 or len(X) != n:
        raise ValueError("X and y must have the same number of rows")
    if not (0.0 < tau < 1.0):
        raise ValueError("tau must be in (0, 1)")
    design = [[1.0] + list(map(float, row)) for row in X] if add_intercept \
        else [list(map(float, row)) for row in X]
    p = len(design[0])
    if n <= p:
        raise ValueError("need more observations than parameters")

    # Start from the equal-weight (OLS) fit.
    beta = _wls(design, y, [1.0] * n, p, n)
    eps = 1e-6
    for _ in range(max_iter):
        weights = []
        for i in range(n):
            fitted = sum(design[i][k] * beta[k] for k in range(p))
            r = y[i] - fitted
            # Asymmetric absolute-loss weight: tau above the fit, 1-tau below.
            asym = tau if r >= 0.0 else (1.0 - tau)
            weights.append(asym / max(abs(r), eps))
        new_beta = _wls(design, y, weights, p, n)
        if max(abs(new_beta[k] - beta[k]) for k in range(p)) < tol:
            beta = new_beta
            break
        beta = new_beta
    return beta


def _wls(design, y, weights, p, n):
    """Weighted least squares: solve (X' W X) beta = X' W y."""
    xtx = [[sum(weights[t] * design[t][i] * design[t][j] for t in range(n))
            for j in range(p)] for i in range(p)]
    xty = [sum(weights[t] * design[t][i] * y[t] for t in range(n)) for i in range(p)]
    inv = _invert(xtx)
    return [sum(inv[i][j] * xty[j] for j in range(p)) for i in range(p)]
