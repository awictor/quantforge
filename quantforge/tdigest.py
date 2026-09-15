"""t-digest: mergeable streaming quantiles with high accuracy in the tails.

The P2 estimator (:mod:`quantforge.streaming_quantile`) tracks one quantile; DDSketch
(:mod:`quantforge.ddsketch`) gives relative-error quantiles on positive data. The t-digest
(Dunning & Ertl) is the workhorse for *arbitrary* quantiles over unbounded streams: it
summarizes the data by a small set of centroids whose sizes are bounded by a scale function
``k(q)`` that keeps clusters tiny near ``q = 0`` and ``q = 1``, so extreme quantiles (p99,
p999) stay accurate while the middle is compressed. Digests merge exactly, which makes them the
standard sketch for distributed percentile aggregation.

This implements the merging variant: buffer incoming points, then sort-and-merge into centroids
under the ``k1`` scale-function size limit. Pure standard library.
"""

import math


class TDigest:
    """A t-digest accumulating a stream of values for approximate quantile queries.

    ``compression`` (delta) trades accuracy for size: larger keeps more centroids and is more
    accurate. Feed values with :meth:`add` (or :meth:`add_all`), read quantiles with
    :meth:`quantile` / :meth:`cdf`, and combine digests with :meth:`merge`.
    """

    def __init__(self, compression=100.0):
        self.compression = float(compression)
        self._centroids = []      # list of [mean, weight], kept sorted by mean
        self._count = 0.0
        self._buffer = []

    def add(self, x, weight=1.0):
        """Add a single value (optionally weighted)."""
        self._buffer.append((float(x), float(weight)))
        if len(self._buffer) > 1000:
            self._flush()
        return self

    def add_all(self, xs):
        for x in xs:
            self.add(x)
        return self

    def _flush(self):
        if not self._buffer:
            return
        # combine existing centroids with the buffer, sort, and re-cluster
        pts = [(c[0], c[1]) for c in self._centroids] + self._buffer
        self._buffer = []
        pts.sort(key=lambda p: p[0])
        total = sum(w for _, w in pts)
        self._count = total
        merged = []
        q0 = 0.0
        cur_mean, cur_w = pts[0]
        for mean, w in pts[1:]:
            q_limit = self._q_limit(q0, total)
            if cur_w + w <= q_limit * total:
                # absorb into current centroid (weighted mean)
                cur_mean = (cur_mean * cur_w + mean * w) / (cur_w + w)
                cur_w += w
            else:
                merged.append([cur_mean, cur_w])
                q0 += cur_w / total
                cur_mean, cur_w = mean, w
        merged.append([cur_mean, cur_w])
        self._centroids = merged

    def _q_limit(self, q0, total):
        # size bound from the k1 scale function: solve k(q0)+1 -> q_limit
        k0 = self._k(q0)
        return self._k_inv(k0 + 1.0) - q0

    def _k(self, q):
        # k1 scale function
        return self.compression / (2.0 * math.pi) * math.asin(2.0 * q - 1.0)

    def _k_inv(self, k):
        return (math.sin(k * 2.0 * math.pi / self.compression) + 1.0) / 2.0

    def merge(self, other):
        """Merge another :class:`TDigest` into this one (exact, order-independent)."""
        self._flush()
        other._flush()
        self._buffer = [(c[0], c[1]) for c in other._centroids]
        self._flush()
        return self

    @property
    def count(self):
        self._flush()
        return self._count

    def quantile(self, q):
        """Approximate value at quantile ``q`` in ``[0, 1]`` (linear interpolation on centroids)."""
        self._flush()
        if not self._centroids:
            raise ValueError("empty digest")
        if q <= 0:
            return self._centroids[0][0]
        if q >= 1:
            return self._centroids[-1][0]
        target = q * self._count
        cum = 0.0
        for i, (mean, w) in enumerate(self._centroids):
            # centroid i covers cumulative weight [cum, cum + w], centered at cum + w/2
            center = cum + w / 2.0
            if target <= center:
                if i == 0:
                    return mean
                pm, pw = self._centroids[i - 1]
                pcenter = (cum - pw) + pw / 2.0
                frac = (target - pcenter) / (center - pcenter)
                return pm + frac * (mean - pm)
            cum += w
        return self._centroids[-1][0]

    def cdf(self, x):
        """Approximate cumulative probability ``P(X <= x)``."""
        self._flush()
        if not self._centroids:
            raise ValueError("empty digest")
        cum = 0.0
        for i, (mean, w) in enumerate(self._centroids):
            center = cum + w / 2.0
            if x < mean:
                if i == 0:
                    return 0.0
                pm, pw = self._centroids[i - 1]
                pcenter = (cum - pw) + pw / 2.0
                frac = (x - pm) / (mean - pm) if mean != pm else 0.0
                return (pcenter + frac * (center - pcenter)) / self._count
            cum += w
        return 1.0
