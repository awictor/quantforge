"""Sparse table: O(1) idempotent range queries on a static array.

For an array that never changes, a sparse table precomputes the answer for every
power-of-two-length block in ``O(n log n)`` time and space, then answers a range query in
``O(1)`` by combining two overlapping blocks that cover the range. This overlap is only
valid for *idempotent* combiners -- ones where ``f(x, x) = x`` and double-counting is
harmless -- so it fits minimum, maximum, and gcd (but not sum; use a Fenwick/segment tree
for that). Pure standard library.
"""

import math


class SparseTable:
    """Static array supporting ``O(1)`` idempotent range queries (min/max/gcd/...).

    Build from a sequence and an associative *idempotent* ``combine`` (default `min`).
    ``query(lo, hi)`` returns the combined value over the inclusive index range
    ``[lo, hi]`` in constant time. The array is fixed after construction.
    """

    __slots__ = ("_combine", "_table", "_n", "_log")

    def __init__(self, values, combine=min):
        data = list(values)
        self._n = len(data)
        self._combine = combine
        if self._n == 0:
            self._table = []
            self._log = [0]
            return
        # log[i] = floor(log2(i)) for block-length lookup
        self._log = [0] * (self._n + 1)
        for i in range(2, self._n + 1):
            self._log[i] = self._log[i // 2] + 1
        max_level = self._log[self._n] + 1
        # table[k][i] = combine of the block starting at i of length 2^k
        table = [data[:]]
        for k in range(1, max_level):
            span = 1 << k
            half = 1 << (k - 1)
            prev = table[k - 1]
            row = []
            for i in range(self._n - span + 1):
                row.append(combine(prev[i], prev[i + half]))
            table.append(row)
        self._table = table

    def __len__(self):
        return self._n

    def query(self, lo, hi):
        """Combine over the inclusive range ``[lo, hi]`` in ``O(1)``."""
        if not 0 <= lo <= hi < self._n:
            raise IndexError("range out of bounds")
        k = self._log[hi - lo + 1]
        span = 1 << k
        left = self._table[k][lo]
        right = self._table[k][hi - span + 1]
        return self._combine(left, right)


def range_min_query(values):
    """Convenience: a `SparseTable` answering range-*minimum* queries."""
    return SparseTable(values, min)


def range_max_query(values):
    """Convenience: a `SparseTable` answering range-*maximum* queries."""
    return SparseTable(values, max)


def range_gcd_query(values):
    """Convenience: a `SparseTable` answering range-*gcd* queries."""
    return SparseTable(values, math.gcd)
