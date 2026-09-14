"""Lazy-propagation segment tree: range updates and range queries in O(log n).

The plain segment tree does point updates; adding *lazy propagation* lets a whole range be
updated at once by tagging interior nodes and pushing the tag down only when a query
descends through them. This supports "add a value to every element of [l, r)" together with
"sum / min / max over [l, r)", each in ``O(log n)``. Three ready-made variants are provided
(sum, min, max under range-add); the class is reusable with a custom aggregate. Pure
standard library.
"""


class LazySegmentTree:
    """Segment tree with range-add updates and a range aggregate query.

    ``mode`` is ``"sum"``, ``"min"``, or ``"max"``. Build from an initial list; then
    ``update(lo, hi, delta)`` adds ``delta`` to every index in ``[lo, hi)`` and
    ``query(lo, hi)`` returns the aggregate over ``[lo, hi)`` -- both ``O(log n)``.
    """

    __slots__ = ("_n", "_mode", "_size", "_tree", "_lazy", "_ident")

    def __init__(self, values, mode="sum"):
        if mode not in ("sum", "min", "max"):
            raise ValueError("mode must be 'sum', 'min', or 'max'")
        data = list(values)
        self._n = len(data)
        self._mode = mode
        self._ident = 0 if mode == "sum" else (float("inf") if mode == "min" else float("-inf"))
        size = 1
        while size < max(1, self._n):
            size <<= 1
        self._size = size
        self._tree = [self._ident] * (2 * size)
        self._lazy = [0] * (2 * size)
        for i, v in enumerate(data):
            self._tree[size + i] = v
        for i in range(size - 1, 0, -1):
            self._tree[i] = self._combine(self._tree[2 * i], self._tree[2 * i + 1])

    def _combine(self, a, b):
        if self._mode == "sum":
            return a + b
        if self._mode == "min":
            return a if a < b else b
        return a if a > b else b

    def _apply(self, node, delta, length):
        self._lazy[node] += delta
        if self._mode == "sum":
            self._tree[node] += delta * length
        else:
            self._tree[node] += delta

    def _push_down(self, node, length):
        if self._lazy[node]:
            half = length // 2
            self._apply(2 * node, self._lazy[node], half)
            self._apply(2 * node + 1, self._lazy[node], half)
            self._lazy[node] = 0

    def __len__(self):
        return self._n

    def update(self, lo, hi, delta):
        """Add ``delta`` to every index in ``[lo, hi)`` (half-open)."""
        if not (0 <= lo <= hi <= self._n):
            raise IndexError("range out of bounds")
        if lo == hi:
            return
        self._update(1, 0, self._size, lo, hi, delta)

    def _update(self, node, nlo, nhi, lo, hi, delta):
        if lo <= nlo and nhi <= hi:
            self._apply(node, delta, nhi - nlo)
            return
        self._push_down(node, nhi - nlo)
        mid = (nlo + nhi) // 2
        if lo < mid:
            self._update(2 * node, nlo, mid, lo, hi, delta)
        if hi > mid:
            self._update(2 * node + 1, mid, nhi, lo, hi, delta)
        self._tree[node] = self._combine(self._tree[2 * node], self._tree[2 * node + 1])

    def query(self, lo, hi):
        """Aggregate over ``[lo, hi)`` (half-open). Empty range returns the identity."""
        if not (0 <= lo <= hi <= self._n):
            raise IndexError("range out of bounds")
        if lo == hi:
            return self._ident
        return self._query(1, 0, self._size, lo, hi)

    def _query(self, node, nlo, nhi, lo, hi):
        if lo <= nlo and nhi <= hi:
            return self._tree[node]
        self._push_down(node, nhi - nlo)
        mid = (nlo + nhi) // 2
        if hi <= mid:
            return self._query(2 * node, nlo, mid, lo, hi)
        if lo >= mid:
            return self._query(2 * node + 1, mid, nhi, lo, hi)
        return self._combine(
            self._query(2 * node, nlo, mid, lo, mid),
            self._query(2 * node + 1, mid, nhi, mid, hi),
        )
