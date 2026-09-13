"""Nonparametric two-sample scale tests: Ansari-Bradley and Mood.

Rank-sum tests (Mann-Whitney, Wilcoxon) detect a *location* shift but are blind to a
pure difference in *spread*. These two test equal dispersion without assuming
normality, by ranking the pooled sample from the outside in and summing the scores of
one group:

  * ``ansari_bradley_test`` -- assigns each pooled observation the score
    ``min(rank, N + 1 - rank)`` (extremes get 1, the center gets the largest score),
    so a more concentrated sample earns higher scores. Assumes the two medians are
    (about) equal.
  * ``mood_test`` -- scores each observation by its squared deviation from the average
    pooled rank ``(rank - (N+1)/2)^2``, summing over one group; a more dispersed
    sample pushes its observations to the rank extremes and inflates the sum.

Both return a normal-approximation two-sided p-value. Univariate; pure standard
library.
"""

import math


def _midranks(pool_sorted, sample):
    """Midranks (1-based, ties averaged) of each value of ``sample`` in the pool."""
    n = len(pool_sorted)
    out = []
    for v in sample:
        lo, hi = 0, n
        while lo < hi:
            mid = (lo + hi) // 2
            if pool_sorted[mid] < v:
                lo = mid + 1
            else:
                hi = mid
        first = lo
        lo2, hi2 = 0, n
        while lo2 < hi2:
            mid = (lo2 + hi2) // 2
            if pool_sorted[mid] <= v:
                lo2 = mid + 1
            else:
                hi2 = mid
        last = lo2 - 1
        out.append((first + 1 + last + 1) / 2.0)
    return out


def ansari_bradley_test(x, y):
    """Ansari-Bradley two-sample test of equal dispersion.

    Scores the pooled ranks from the outside in (``min(r, N+1-r)``) and sums the
    scores over ``x``. Returns a dict with the ``statistic`` (that sum), the ``z``
    normal approximation and the two-sided ``p_value``. Assumes the two samples share
    a location; a smaller ``x`` spread pushes ``x`` toward the center (higher scores).
    """
    m = len(x)
    n = len(y)
    if m == 0 or n == 0:
        raise ValueError("both samples must be non-empty")
    big_n = m + n
    pool = sorted(list(x) + list(y))
    rx = _midranks(pool, x)
    scores = [min(r, big_n + 1 - r) for r in rx]
    stat = sum(scores)

    # Null mean/variance of the Ansari-Bradley statistic (exact for the no-tie case).
    if big_n % 2 == 0:
        mean = m * (big_n + 2) / 4.0
        var = m * n * (big_n + 2) * (big_n - 2) / (48.0 * (big_n - 1))
    else:
        mean = m * (big_n + 1) ** 2 / (4.0 * big_n)
        var = m * n * (big_n + 1) * (3 + big_n * big_n) / (48.0 * big_n * big_n)
    if var <= 0:
        z = 0.0
    else:
        z = (stat - mean) / math.sqrt(var)
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return {"statistic": stat, "z": z, "p_value": p}


def mood_test(x, y):
    """Mood two-sample test of equal dispersion.

    Sums the squared rank deviations ``(r - (N+1)/2)^2`` over ``x``. Returns a dict
    with the ``statistic``, the ``z`` normal approximation and the two-sided
    ``p_value``. A more dispersed ``x`` sends its values to the rank extremes and
    raises the statistic. Assumes a common location.
    """
    m = len(x)
    n = len(y)
    if m == 0 or n == 0:
        raise ValueError("both samples must be non-empty")
    big_n = m + n
    pool = sorted(list(x) + list(y))
    rx = _midranks(pool, x)
    centre = (big_n + 1) / 2.0
    stat = sum((r - centre) ** 2 for r in rx)

    mean = m * (big_n * big_n - 1) / 12.0
    var = m * n * (big_n + 1) * (big_n * big_n - 4) / 180.0
    if var <= 0:
        z = 0.0
    else:
        z = (stat - mean) / math.sqrt(var)
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return {"statistic": stat, "z": z, "p_value": p}
