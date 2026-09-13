"""Distance correlation (Szekely-Rizzo): a dependence measure that sees nonlinearity.

Pearson correlation is blind to nonlinear dependence -- a point on a circle or a
parabola has near-zero Pearson r yet is strongly dependent. Distance correlation
fixes this: it is zero *if and only if* the two variables are independent, for any
kind of relationship. Build the pairwise-distance matrix of each sample, double-center
it (subtract row means, column means, add back the grand mean), and the distance
covariance is the mean product of the two centered matrices:

    dCov^2(X, Y) = mean_{i,j} A_ij B_ij,
    dCor(X, Y)   = dCov(X, Y) / sqrt(dVar(X) dVar(Y)).

``dCor`` lies in ``[0, 1]``: 0 exactly under independence, 1 for a tight linear
relation. Univariate samples here; ``O(n^2)``. Pure standard library.
"""

import math


def _double_centered(sample):
    """Double-centered pairwise |distance| matrix of a 1-D sample."""
    n = len(sample)
    dist = [[abs(sample[i] - sample[j]) for j in range(n)] for i in range(n)]
    row_mean = [sum(row) / n for row in dist]
    grand = sum(row_mean) / n
    col_mean = row_mean  # symmetric matrix -> column means equal row means
    return [[dist[i][j] - row_mean[i] - col_mean[j] + grand
             for j in range(n)] for i in range(n)], n


def distance_covariance(x, y):
    """Distance covariance ``dCov(x, y)`` (the square root of the mean product).

    Zero if and only if ``x`` and ``y`` are independent (in the population). Returns a
    non-negative value on the same scale as the data.
    """
    if len(x) != len(y):
        raise ValueError("x and y must have equal length")
    if len(x) < 2:
        raise ValueError("need at least 2 points")
    A, n = _double_centered(x)
    B, _ = _double_centered(y)
    s = sum(A[i][j] * B[i][j] for i in range(n) for j in range(n)) / (n * n)
    return math.sqrt(s) if s > 0 else 0.0


def distance_variance(x):
    """Distance variance ``dVar(x) = dCov(x, x)``; zero only for a constant sample."""
    A, n = _double_centered(x)
    s = sum(A[i][j] * A[i][j] for i in range(n) for j in range(n)) / (n * n)
    return math.sqrt(s) if s > 0 else 0.0


def distance_correlation(x, y):
    """Distance correlation ``dCor(x, y)`` in ``[0, 1]``.

    ``dCov(x, y) / sqrt(dVar(x) dVar(y))``. Zero exactly under independence (unlike
    Pearson, this holds for nonlinear dependence too), 1 for a tight linear relation.
    Returns 0 if either sample is constant (distance variance zero).
    """
    if len(x) != len(y):
        raise ValueError("x and y must have equal length")
    if len(x) < 2:
        raise ValueError("need at least 2 points")
    A, n = _double_centered(x)
    B, _ = _double_centered(y)
    dcov2 = sum(A[i][j] * B[i][j] for i in range(n) for j in range(n)) / (n * n)
    dvarx = sum(A[i][j] * A[i][j] for i in range(n) for j in range(n)) / (n * n)
    dvary = sum(B[i][j] * B[i][j] for i in range(n) for j in range(n)) / (n * n)
    denom = math.sqrt(dvarx * dvary)
    if denom <= 0:
        return 0.0
    return math.sqrt(max(dcov2, 0.0) / denom)
