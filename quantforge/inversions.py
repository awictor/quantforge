"""Counting inversions and ranking distances.

An *inversion* is a pair of positions out of order -- ``i < j`` but ``a[i] > a[j]``. The
number of inversions measures how far a sequence is from sorted (0 when sorted, ``n(n-1)/2``
when reversed), and a merge sort counts them all in ``O(n log n)`` as a byproduct of the
merge. The same idea gives the Kendall-tau distance between two rankings -- the count of
pairs they order oppositely, i.e. the number of adjacent swaps a bubble sort would need.
Pure standard library.
"""


def count_inversions(values):
    """Number of inversions in ``values`` (pairs ``i < j`` with ``values[i] > values[j]``).

    Uses a merge sort, so ``O(n log n)`` rather than the ``O(n^2)`` brute count. A strictly
    increasing sequence has 0; a strictly decreasing one has ``n(n-1)/2``.
    """
    def sort_count(arr):
        n = len(arr)
        if n <= 1:
            return arr, 0
        mid = n // 2
        left, cl = sort_count(arr[:mid])
        right, cr = sort_count(arr[mid:])
        merged = []
        i = j = 0
        inv = cl + cr
        nl = len(left)
        while i < nl and j < len(right):
            if left[i] <= right[j]:
                merged.append(left[i])
                i += 1
            else:
                merged.append(right[j])
                j += 1
                inv += nl - i          # left[i:] all exceed right[j]
        merged.extend(left[i:])
        merged.extend(right[j:])
        return merged, inv

    _, total = sort_count(list(values))
    return total


def kendall_tau_distance(rank_a, rank_b):
    """Kendall-tau distance: number of pairs ordered oppositely in two rankings.

    ``rank_a`` and ``rank_b`` are sequences of the same items (equal length, same set). The
    distance is the count of pairs ``(x, y)`` whose relative order differs between the two,
    i.e. the inversions of ``rank_b`` reindexed by ``rank_a``'s positions.
    """
    if len(rank_a) != len(rank_b):
        raise ValueError("rankings must have equal length")
    if set(rank_a) != set(rank_b):
        raise ValueError("rankings must contain the same items")
    pos = {item: i for i, item in enumerate(rank_a)}
    if len(pos) != len(rank_a):
        raise ValueError("rankings must not contain duplicates")
    mapped = [pos[item] for item in rank_b]
    return count_inversions(mapped)


def is_sorted(values, strict=False):
    """True if ``values`` is non-decreasing (or strictly increasing when ``strict``)."""
    if strict:
        return all(values[i] < values[i + 1] for i in range(len(values) - 1))
    return all(values[i] <= values[i + 1] for i in range(len(values) - 1))
