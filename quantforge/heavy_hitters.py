"""Streaming frequent-items (heavy hitters): Misra-Gries and Space-Saving.

Finding the items that make up a large share of a stream, using far less memory than a
full frequency table:

  * ``MisraGries`` -- the classic bounded-counter summary. With ``k`` counters, any item
    occurring more than ``n / (k + 1)`` times is guaranteed to be kept, and each stored
    count undercounts the true frequency by at most ``n / (k + 1)``.
  * ``SpaceSaving`` -- Metwally-Agarwal-Abbadi, which keeps ``k`` (item, count, error)
    triples and, on overflow, evicts the current minimum and reuses its slot. Its stored
    count *overcounts* by at most the recorded error, so ``count - error`` is a guaranteed
    lower bound; it tends to rank the true top-k more accurately than Misra-Gries.

Both run in one pass and ``O(k)`` memory. Pure standard library.
"""


class MisraGries:
    """Misra-Gries frequent-items summary with ``k`` counters.

    ``add(item)`` folds one occurrence; ``counts()`` returns the surviving item -> count
    map. Any item with true frequency greater than ``n / (k + 1)`` is guaranteed present,
    and stored counts never exceed the true count (they may undercount by up to the number
    of decrement rounds).
    """

    __slots__ = ("k", "n", "_counts")

    def __init__(self, k):
        if k < 1:
            raise ValueError("k must be at least 1")
        self.k = int(k)
        self.n = 0
        self._counts = {}

    def add(self, item, weight=1):
        """Fold ``weight`` occurrences of ``item`` (weight must be a positive integer)."""
        if weight < 1:
            raise ValueError("weight must be a positive integer")
        self.n += weight
        c = self._counts
        if item in c:
            c[item] += weight
        elif len(c) < self.k:
            c[item] = weight
        else:
            # decrement all counters; drop any that hit zero. Charge the whole weight.
            dec = weight
            # the incoming item effectively cancels against the decrement rounds
            for key in list(c.keys()):
                c[key] -= dec
                if c[key] <= 0:
                    del c[key]
        return self

    def counts(self):
        """Return a copy of the surviving item -> approximate-count map."""
        return dict(self._counts)

    def heavy_hitters(self, threshold):
        """Items whose stored count meets ``threshold`` (a fraction in (0, 1] of ``n``)."""
        if not 0.0 < threshold <= 1.0:
            raise ValueError("threshold must be in (0, 1]")
        cutoff = threshold * self.n
        return {k: v for k, v in self._counts.items() if v >= cutoff}


class SpaceSaving:
    """Space-Saving top-k summary (Metwally-Agarwal-Abbadi) with ``k`` counters.

    ``add(item)`` folds one occurrence; ``top(m)`` returns the ``m`` highest
    ``(item, count)`` pairs. Each slot tracks an ``error`` -- the maximum it may overcount
    -- so ``count - error`` is a guaranteed lower bound on the true frequency. When a new
    item arrives and every slot is used, the slot with the smallest count is evicted and
    its count becomes the new item's starting error.
    """

    __slots__ = ("k", "n", "_counts", "_error")

    def __init__(self, k):
        if k < 1:
            raise ValueError("k must be at least 1")
        self.k = int(k)
        self.n = 0
        self._counts = {}
        self._error = {}

    def add(self, item, weight=1):
        """Fold ``weight`` occurrences of ``item`` (weight must be a positive integer)."""
        if weight < 1:
            raise ValueError("weight must be a positive integer")
        self.n += weight
        c = self._counts
        if item in c:
            c[item] += weight
        elif len(c) < self.k:
            c[item] = weight
            self._error[item] = 0
        else:
            # evict the current minimum-count slot and reuse it for the new item
            victim = min(c, key=c.__getitem__)
            min_count = c[victim]
            del c[victim]
            del self._error[victim]
            c[item] = min_count + weight
            self._error[item] = min_count
        return self

    def counts(self):
        """Return a copy of the tracked item -> count map."""
        return dict(self._counts)

    def errors(self):
        """Return a copy of the item -> maximum-overcount map."""
        return dict(self._error)

    def top(self, m):
        """Return the ``m`` highest ``(item, count)`` pairs, descending by count."""
        if m < 0:
            raise ValueError("m must be non-negative")
        ranked = sorted(self._counts.items(), key=lambda kv: kv[1], reverse=True)
        return ranked[:m]

    def guaranteed(self, m):
        """Top ``m`` as ``(item, lower_bound)`` where ``lower_bound = count - error``."""
        ranked = self.top(m)
        return [(it, cnt - self._error[it]) for it, cnt in ranked]
