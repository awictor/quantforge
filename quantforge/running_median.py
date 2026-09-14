"""Exact running median of a stream via two heaps.

`P2Quantile` estimates a quantile in constant memory; when the *exact* median of everything
seen so far is needed, the two-heap method keeps a max-heap of the lower half and a min-heap
of the upper half, rebalanced so their sizes differ by at most one. Each insert is
``O(log n)`` and the median is the top of the larger heap (or the average of both tops when
the counts are equal). Memory grows with the stream, unlike P2. Pure standard library.
"""

import heapq


class RunningMedian:
    """Maintains the exact median of all values inserted so far.

    ``push(x)`` adds a value in ``O(log n)``; ``median()`` returns the current median in
    ``O(1)``. ``lower`` is a max-heap (stored as negated values) of the smaller half; ``upper``
    a min-heap of the larger half.
    """

    __slots__ = ("_lower", "_upper")

    def __init__(self, values=None):
        self._lower = []      # max-heap via negation: the smaller half
        self._upper = []      # min-heap: the larger half
        if values is not None:
            for v in values:
                self.push(v)

    def __len__(self):
        return len(self._lower) + len(self._upper)

    def push(self, x):
        """Insert ``x`` and keep the two heaps balanced."""
        x = float(x)
        if not self._lower or x <= -self._lower[0]:
            heapq.heappush(self._lower, -x)
        else:
            heapq.heappush(self._upper, x)
        # rebalance so |lower| - |upper| in {0, 1}
        if len(self._lower) > len(self._upper) + 1:
            heapq.heappush(self._upper, -heapq.heappop(self._lower))
        elif len(self._upper) > len(self._lower):
            heapq.heappush(self._lower, -heapq.heappop(self._upper))
        return self

    def median(self):
        """Current median; raises if empty. Averages the two middle values on even counts."""
        n = len(self)
        if n == 0:
            raise ValueError("median of an empty stream")
        if len(self._lower) > len(self._upper):
            return -self._lower[0]
        return (-self._lower[0] + self._upper[0]) / 2.0

    def count(self):
        """Number of values seen."""
        return len(self)
