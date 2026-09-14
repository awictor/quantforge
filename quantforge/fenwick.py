"""Fenwick (binary indexed) tree and a generic segment tree.

Data structures for fast prefix/range queries with updates. The Fenwick tree gives
``O(log n)`` point updates and prefix sums (hence range sums) in tight space -- the go-to
for cumulative-frequency and running-total problems. The segment tree generalizes to any
associative combine (sum, min, max) with ``O(log n)`` point update and range query. Pure
standard library.
"""


class FenwickTree:
    """Binary indexed tree for point updates and prefix/range sums (``O(log n)`` each).

    Construct from a size (all zeros) or an initial list. ``update(i, delta)`` adds to
    element ``i``; ``prefix_sum(i)`` is the sum of elements ``0..i``; ``range_sum(lo, hi)``
    the inclusive sum over ``[lo, hi]``. Zero-indexed.
    """

    def __init__(self, size_or_values):
        if isinstance(size_or_values, int):
            self.n = size_or_values
            self.tree = [0.0] * (self.n + 1)
        else:
            values = list(size_or_values)
            self.n = len(values)
            self.tree = [0.0] * (self.n + 1)
            for i, v in enumerate(values):
                self.update(i, v)

    def update(self, i, delta):
        """Add ``delta`` to element ``i`` (0-indexed)."""
        if not (0 <= i < self.n):
            raise IndexError("index out of range")
        i += 1
        while i <= self.n:
            self.tree[i] += delta
            i += i & (-i)

    def prefix_sum(self, i):
        """Sum of elements ``0..i`` inclusive; ``prefix_sum(-1) == 0``."""
        if i < -1 or i >= self.n:
            raise IndexError("index out of range")
        s = 0.0
        i += 1
        while i > 0:
            s += self.tree[i]
            i -= i & (-i)
        return s

    def range_sum(self, lo, hi):
        """Sum over the inclusive range ``[lo, hi]``."""
        if lo > hi:
            raise ValueError("require lo <= hi")
        return self.prefix_sum(hi) - self.prefix_sum(lo - 1)


class SegmentTree:
    """Segment tree for range queries under an associative ``combine`` with point updates.

    ``combine`` defaults to ``+`` (range sum); pass ``min``/``max`` (with the matching
    ``identity``) for range-minimum/maximum. ``update(i, value)`` sets element ``i``;
    ``query(lo, hi)`` folds ``combine`` over ``[lo, hi]`` inclusive. ``O(log n)`` each.
    """

    def __init__(self, values, combine=None, identity=0.0):
        self.combine = combine if combine is not None else (lambda a, b: a + b)
        self.identity = identity
        data = list(values)
        self.n = len(data)
        if self.n == 0:
            raise ValueError("need at least one value")
        self.size = 1
        while self.size < self.n:
            self.size *= 2
        self.tree = [identity] * (2 * self.size)
        for i, v in enumerate(data):
            self.tree[self.size + i] = v
        for i in range(self.size - 1, 0, -1):
            self.tree[i] = self.combine(self.tree[2 * i], self.tree[2 * i + 1])

    def update(self, i, value):
        """Set element ``i`` (0-indexed) to ``value`` and repair the tree."""
        if not (0 <= i < self.n):
            raise IndexError("index out of range")
        pos = self.size + i
        self.tree[pos] = value
        pos //= 2
        while pos >= 1:
            self.tree[pos] = self.combine(self.tree[2 * pos], self.tree[2 * pos + 1])
            pos //= 2

    def query(self, lo, hi):
        """Fold ``combine`` over the inclusive range ``[lo, hi]``."""
        if lo > hi or lo < 0 or hi >= self.n:
            raise ValueError("invalid query range")
        res = self.identity
        lo += self.size
        hi += self.size + 1
        while lo < hi:
            if lo & 1:
                res = self.combine(res, self.tree[lo])
                lo += 1
            if hi & 1:
                hi -= 1
                res = self.combine(res, self.tree[hi])
            lo //= 2
            hi //= 2
        return res
