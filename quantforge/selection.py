"""Selection: k-th order statistic and top-k without a full sort.

Finding the k-th smallest element (or the top-k) does not need an ``O(n log n)`` sort.
Quickselect with the median-of-medians pivot runs in guaranteed linear time even on
adversarial input. Provides the k-th order statistic, the median, and the k smallest or
largest elements. Pure standard library.
"""


def _median_of_medians(a):
    """A pivot value guaranteed to be near the median (median-of-medians)."""
    if len(a) <= 5:
        return sorted(a)[len(a) // 2]
    medians = [sorted(a[i:i + 5])[len(a[i:i + 5]) // 2] for i in range(0, len(a), 5)]
    return _median_of_medians(medians)


def kth_smallest(values, k):
    """The ``k``-th smallest element (0-indexed) via linear-time quickselect.

    ``k`` in ``[0, len(values) - 1]``; ``k=0`` is the minimum. Uses the median-of-medians
    pivot for a guaranteed ``O(n)`` worst case. Does not modify the input.
    """
    a = list(values)
    n = len(a)
    if n == 0:
        raise ValueError("values must be non-empty")
    if not (0 <= k < n):
        raise ValueError("k out of range")
    while True:
        if len(a) == 1:
            return a[0]
        pivot = _median_of_medians(a)
        lows = [x for x in a if x < pivot]
        highs = [x for x in a if x > pivot]
        pivots = [x for x in a if x == pivot]
        if k < len(lows):
            a = lows
        elif k < len(lows) + len(pivots):
            return pivot
        else:
            k -= len(lows) + len(pivots)
            a = highs


def median(values):
    """Median of ``values`` (average of the two middle elements for even length)."""
    n = len(values)
    if n == 0:
        raise ValueError("values must be non-empty")
    if n % 2 == 1:
        return kth_smallest(values, n // 2)
    lo = kth_smallest(values, n // 2 - 1)
    hi = kth_smallest(values, n // 2)
    return 0.5 * (lo + hi)


def top_k(values, k, largest=True):
    """The ``k`` largest (or smallest) elements of ``values``, sorted.

    ``largest=True`` returns the top ``k`` in descending order; ``largest=False`` the
    bottom ``k`` ascending. ``k`` is clamped to the list length. Uses quickselect to find
    the threshold, then sorts only the ``k`` selected elements.
    """
    a = list(values)
    n = len(a)
    if k < 0:
        raise ValueError("k must be non-negative")
    if k == 0:
        return []
    if k >= n:
        return sorted(a, reverse=largest)
    # Index of the threshold order statistic.
    idx = n - k if largest else k - 1
    threshold = kth_smallest(a, idx)
    if largest:
        chosen = [x for x in a if x > threshold]
        chosen += [threshold] * (k - len(chosen))
        return sorted(chosen, reverse=True)
    else:
        chosen = [x for x in a if x < threshold]
        chosen += [threshold] * (k - len(chosen))
        return sorted(chosen)
