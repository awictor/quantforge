"""Theil-Sen robust linear regression.

Ordinary least squares is pulled hard by outliers -- a single bad point tips the
fitted line. The Theil-Sen estimator instead takes the *median* of the slopes of
all pairs of points:

    slope = median_{i < j} (y_j - y_i) / (x_j - x_i),
    intercept = median_i (y_i - slope * x_i).

The median gives it a breakdown point of about 29.3% -- up to that fraction of the
data can be arbitrarily corrupted before the estimate blows up -- while remaining
exact on clean linear data and nearly as efficient as OLS on Gaussian noise. It is
the standard robust slope for trend detection. Pure standard library.
"""


def _median(values):
    """Median of a non-empty list (average of the two middle order statistics)."""
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2 == 1:
        return s[mid]
    return 0.5 * (s[mid - 1] + s[mid])


def theil_sen(x, y):
    """Theil-Sen robust regression slope and intercept.

    Parameters
    ----------
    x, y : sequence of float
        Paired observations of equal length (at least 2). Pairs sharing an ``x``
        value are skipped (their slope is undefined).

    Returns
    -------
    (slope, intercept) : (float, float)
        The median pairwise slope and the median residual intercept. On exactly
        collinear data this reproduces the generating line; under heavy-tailed
        contamination it stays close to the clean fit where OLS is dragged away.
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    n = len(x)
    if n < 2:
        raise ValueError("need at least 2 points")

    slopes = []
    for i in range(n):
        for j in range(i + 1, n):
            dx = x[j] - x[i]
            if dx != 0.0:
                slopes.append((y[j] - y[i]) / dx)
    if not slopes:
        raise ValueError("all x values are identical; slope undefined")
    slope = _median(slopes)
    intercept = _median([y[i] - slope * x[i] for i in range(n)])
    return slope, intercept
