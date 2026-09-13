"""Deming regression: errors-in-variables (total least squares) line fitting.

Ordinary least squares assumes ``x`` is measured without error and minimizes only
vertical residuals; when *both* variables carry measurement error (two instruments,
two assays) that biases the slope toward zero. Deming regression accounts for error in
both, minimizing the sum of squared distances weighted by the error-variance ratio
``lambda = var(error_x) / var(error_y)``. The closed-form slope is

    b = (Syy - lambda Sxx + sqrt((Syy - lambda Sxx)^2 + 4 lambda Sxy^2)) / (2 Sxy),

with ``Sxx, Syy, Sxy`` the (co)variance sums. ``lambda -> infinity`` recovers OLS of
``y`` on ``x``; ``lambda = 1`` is orthogonal / total-least-squares regression. Pure
standard library.
"""

import math


def deming_regression(x, y, lam=1.0):
    """Deming (errors-in-variables) regression of ``y`` on ``x``.

    ``lam`` is the ratio of the error variances ``var(eps_x) / var(eps_y)`` (default 1,
    i.e. orthogonal regression). Returns ``(slope, intercept)`` for the fitted line
    ``y = slope * x + intercept``. Requires a non-zero ``Sxy``.
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have equal length")
    if n < 2:
        raise ValueError("need at least 2 points")
    if lam <= 0:
        raise ValueError("lambda must be positive")

    mx = sum(x) / n
    my = sum(y) / n
    sxx = sum((xi - mx) ** 2 for xi in x) / (n - 1)
    syy = sum((yi - my) ** 2 for yi in y) / (n - 1)
    sxy = sum((x[i] - mx) * (y[i] - my) for i in range(n)) / (n - 1)
    if sxy == 0:
        raise ValueError("Sxy is zero; slope undefined for Deming regression")

    disc = (syy - lam * sxx) ** 2 + 4.0 * lam * sxy * sxy
    slope = (syy - lam * sxx + math.sqrt(disc)) / (2.0 * sxy)
    intercept = my - slope * mx
    return slope, intercept


def orthogonal_regression(x, y):
    """Orthogonal (total-least-squares) regression: Deming with ``lambda = 1``.

    Minimizes the perpendicular distances to the line, treating ``x`` and ``y``
    symmetrically. Returns ``(slope, intercept)``.
    """
    return deming_regression(x, y, lam=1.0)
