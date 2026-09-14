"""Order-statistic tree over a fixed value universe (Fenwick-backed).

A multiset that answers order queries in ``O(log U)``: insert or remove a value, ask its
rank (how many stored values are smaller), select the ``k``-th smallest, or count how many
lie below a threshold or inside a range. Built on a binary indexed tree over the value
universe supplied at construction, so it fits ranking, running percentiles, and inversion
counting on a known key range. Pure standard library.
"""

import bisect


class OrderStatisticTree:
    """Dynamic multiset with rank/select over a fixed sorted value universe.

    Construct with the universe of possible values (deduplicated and sorted internally).
    ``add``/``remove`` adjust multiplicities; ``rank(x)`` counts stored values ``< x``;
    ``select(k)`` returns the ``k``-th smallest (0-indexed); ``count_less``/``count_range``
    answer threshold and interval counts. All queries are ``O(log U)``.
    """

    __slots__ = ("_values", "_index", "_tree", "_n", "_size")

    def __init__(self, universe):
        self._values = sorted(set(universe))
        self._index = {v: i + 1 for i, v in enumerate(self._values)}  # 1-based BIT positions
        self._n = len(self._values)
        self._tree = [0] * (self._n + 1)
        self._size = 0

    def __len__(self):
        return self._size

    def _update(self, pos, delta):
        while pos <= self._n:
            self._tree[pos] += delta
            pos += pos & -pos

    def _prefix(self, pos):
        s = 0
        while pos > 0:
            s += self._tree[pos]
            pos -= pos & -pos
        return s

    def add(self, value, count=1):
        """Insert ``count`` copies of ``value`` (must be in the universe)."""
        if value not in self._index:
            raise ValueError("value not in the universe")
        if count < 1:
            raise ValueError("count must be positive")
        self._update(self._index[value], count)
        self._size += count
        return self

    def remove(self, value, count=1):
        """Remove ``count`` copies of ``value``; raises if fewer are present."""
        if value not in self._index:
            raise ValueError("value not in the universe")
        pos = self._index[value]
        present = self._prefix(pos) - self._prefix(pos - 1)
        if count > present:
            raise ValueError("cannot remove more copies than are present")
        self._update(pos, -count)
        self._size -= count
        return self

    def count_less(self, x):
        """Number of stored values strictly less than ``x``."""
        # rank of x = prefix count up to the largest universe value < x
        i = bisect.bisect_left(self._values, x)
        return self._prefix(i)

    def rank(self, x):
        """Alias for :meth:`count_less`: how many stored values are ``< x``."""
        return self.count_less(x)

    def count_range(self, lo, hi):
        """Number of stored values in the half-open range ``[lo, hi)``."""
        if lo > hi:
            raise ValueError("lo must not exceed hi")
        return self.count_less(hi) - self.count_less(lo)

    def count(self, value):
        """Multiplicity of ``value`` currently stored."""
        if value not in self._index:
            return 0
        pos = self._index[value]
        return self._prefix(pos) - self._prefix(pos - 1)

    def select(self, k):
        """Return the ``k``-th smallest stored value (0-indexed). Raises if out of range."""
        if not 0 <= k < self._size:
            raise IndexError("k out of range")
        target = k + 1
        pos = 0
        acc = 0
        log = self._n.bit_length()
        step = 1 << log
        while step > 0:
            nxt = pos + step
            if nxt <= self._n and acc + self._tree[nxt] < target:
                pos = nxt
                acc += self._tree[nxt]
            step >>= 1
        return self._values[pos]  # pos is 0-based index of the (pos+1)-th universe value
