"""Hodges-Lehmann robust location and shift estimators.

The Hodges-Lehmann estimator is the median of pairwise averages -- a robust location
estimator with a 29% breakdown point (far above the mean's 0%) yet ~96% efficiency at
the normal, and the point estimate paired with the Wilcoxon signed-rank / rank-sum
tests. Two forms:

  * ``hodges_lehmann_location`` -- one sample: median over all Walsh averages
    ``(x_i + x_j) / 2`` for ``i <= j``. The location a symmetric distribution is
    centered on.
  * ``hodges_lehmann_shift`` -- two samples: median over all pairwise differences
    ``y_j - x_i``. The typical shift between the two, the estimate that inverts the
    Mann-Whitney / Wilcoxon rank-sum test.

Pure standard library; ``O(n^2)`` pair enumeration with an exact median.
"""


def _median_sorted(values):
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2:
        return s[mid]
    return 0.5 * (s[mid - 1] + s[mid])


def hodges_lehmann_location(x):
    """One-sample Hodges-Lehmann estimator: median of the Walsh averages.

    Median of ``(x_i + x_j) / 2`` over all ``i <= j`` (including ``i == j``). Robust
    (29% breakdown) and highly efficient at the normal; estimates the center of a
    symmetric distribution.
    """
    n = len(x)
    if n == 0:
        raise ValueError("sample must be non-empty")
    walsh = [(x[i] + x[j]) / 2.0 for i in range(n) for j in range(i, n)]
    return _median_sorted(walsh)


def hodges_lehmann_shift(x, y):
    """Two-sample Hodges-Lehmann shift: median of all pairwise differences.

    Median of ``y_j - x_i`` over every pair. The robust estimate of the location
    shift between the two samples, consistent with the Wilcoxon rank-sum test (the
    shift for which the test would not reject). Positive means ``y`` is shifted above
    ``x``.
    """
    if len(x) == 0 or len(y) == 0:
        raise ValueError("both samples must be non-empty")
    diffs = [yj - xi for xi in x for yj in y]
    return _median_sorted(diffs)
