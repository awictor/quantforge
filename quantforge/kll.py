"""KLL sketch: streaming quantiles with a provable rank-error guarantee.

The t-digest (:mod:`quantforge.tdigest`) is accurate in practice but has no worst-case bound.
The KLL sketch (Karnin, Lang & Liberty 2016) is the theoretically optimal comparison-based
quantile sketch: it answers any rank query within an additive error of about ``epsilon * n``
using ``O((1/epsilon) log log ...)`` space, and it merges. It works by keeping a hierarchy of
*compactors* -- sorted buffers, one per level; when a level fills, it sorts, keeps every other
element (a random even/odd offset), and promotes them to the next level, where each carries
double the weight.

A seeded :class:`quantforge.pcg.PCG32` stream drives the compaction coin flips so runs are
reproducible. Pure standard library.
"""

from .pcg import PCG32


class KLL:
    """A KLL quantile sketch. Larger ``k`` -> smaller rank error (``~ c / k``) and more memory.

    Feed values with :meth:`add`/:meth:`add_all`, query with :meth:`rank`, :meth:`quantile`,
    :meth:`cdf`; combine sketches with :meth:`merge`.
    """

    def __init__(self, k=200, c=2.0 / 3.0, seed=12345):
        self.k = int(k)
        self.c = c
        self._rng = PCG32(seed)
        self._levels = [[]]      # levels[h] is the compactor at height h
        self._n = 0

    def add(self, x):
        self._levels[0].append(float(x))
        self._n += 1
        self._compress()
        return self

    def add_all(self, xs):
        for x in xs:
            self.add(x)
        return self

    def _capacity(self, height, num_levels):
        # capacity shrinks geometrically for lower levels; top level gets full k
        depth = num_levels - height - 1
        cap = int(self.k * (self.c ** depth))
        return max(cap, 2)

    def _compress(self):
        h = 0
        while h < len(self._levels):
            num_levels = len(self._levels)
            cap = self._capacity(h, num_levels)
            if len(self._levels[h]) < cap:
                h += 1
                continue
            # compact level h into level h+1
            if h + 1 == len(self._levels):
                self._levels.append([])
            buf = sorted(self._levels[h])
            offset = 1 if self._rng.random() < 0.5 else 0
            promoted = buf[offset::2]     # keep every other element
            self._levels[h] = []
            self._levels[h + 1].extend(promoted)
            h += 1

    def _weighted_items(self):
        # (value, weight) pairs; weight = 2^height
        items = []
        for height, buf in enumerate(self._levels):
            w = 1 << height
            for v in buf:
                items.append((v, w))
        return items

    def rank(self, x):
        """Estimated number of stream elements ``<= x`` (weighted count over all levels)."""
        return sum(w for v, w in self._weighted_items() if v <= x)

    @property
    def n(self):
        return self._n

    def cdf(self, x):
        """Estimated ``P(X <= x)``."""
        if self._n == 0:
            raise ValueError("empty sketch")
        return self.rank(x) / self._n

    def quantile(self, q):
        """Approximate value at quantile ``q`` (rank ``q * n``), within ~``c/k`` rank error."""
        if self._n == 0:
            raise ValueError("empty sketch")
        items = sorted(self._weighted_items(), key=lambda p: p[0])
        target = q * self._n
        cum = 0.0
        for v, w in items:
            cum += w
            if cum >= target:
                return v
        return items[-1][0]

    def merge(self, other):
        """Merge another :class:`KLL` sketch in (concatenate levels, then re-compress)."""
        while len(self._levels) < len(other._levels):
            self._levels.append([])
        for h, buf in enumerate(other._levels):
            self._levels[h].extend(buf)
        self._n += other._n
        self._compress()
        return self
