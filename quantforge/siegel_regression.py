"""Siegel's repeated-median regression: 50% breakdown robust line fitting.

Theil-Sen takes the median of all pairwise slopes and breaks down at ~29% outliers.
Siegel's (1982) repeated median goes further: for each point ``i`` take the median of
the slopes to every *other* point, then take the median of those per-point medians.
That double median tolerates up to 50% contaminated data -- the highest possible
breakdown for a regression estimator -- while remaining exact on clean linear data. The
intercept is the median of ``y_i - slope * x_i``. Pure standard library; ``O(n^2)``.
"""


def _median(values):
    s = sorted(values)
    n = len(s)
    if n == 0:
        raise ValueError("empty")
    mid = n // 2
    if n % 2:
        return s[mid]
    return 0.5 * (s[mid - 1] + s[mid])


def repeated_median_regression(x, y):
    """Siegel repeated-median regression; returns ``(slope, intercept)``.

    For each point, the median of its slopes to all other points; the overall slope is
    the median of those. 50% breakdown point -- robust to nearly half the data being
    corrupted. Points sharing an ``x`` value contribute no slope for that pair.
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have equal length")
    if n < 2:
        raise ValueError("need at least 2 points")

    per_point = []
    for i in range(n):
        slopes = []
        for j in range(n):
            if j != i and x[j] != x[i]:
                slopes.append((y[j] - y[i]) / (x[j] - x[i]))
        if slopes:
            per_point.append(_median(slopes))
    if not per_point:
        raise ValueError("all x values are identical; slope undefined")
    slope = _median(per_point)
    intercept = _median([y[i] - slope * x[i] for i in range(n)])
    return slope, intercept
