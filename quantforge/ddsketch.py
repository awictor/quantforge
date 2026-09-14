"""DDSketch: a mergeable quantile sketch with a relative-error guarantee.

Unlike a fixed-quantile tracker (`P2Quantile`) or an equiprobable reservoir, DDSketch
(Masson-Rim-Lee 2019) answers *any* quantile after the fact and guarantees a bound on the
*relative* error of the returned value: for accuracy ``alpha``, the reported ``q``-quantile
``x_hat`` satisfies ``|x_hat - x_q| <= alpha * x_q``. It works by mapping each positive
value to a bucket ``i = ceil(log(x) / log(gamma))`` with ``gamma = (1 + alpha)/(1 - alpha)``
and counting bucket populations, so the buckets are logarithmically spaced -- fine where
the data is small, coarse where it is large. Sketches over the same ``alpha`` merge exactly
by adding bucket counts. Positive values only (the usual DDSketch domain: latencies, sizes,
durations). Pure standard library.
"""

import math


class DDSketch:
    """Relative-error quantile sketch over positive values.

    ``add(x)`` folds one positive value; ``quantile(q)`` returns the estimated
    ``q``-quantile with relative error at most ``alpha``. ``merge`` / ``+`` combine two
    sketches built with the same ``alpha``. Memory grows only with the log-ratio of the
    largest to smallest value seen, not with the number of points.
    """

    __slots__ = ("alpha", "_gamma", "_log_gamma", "_buckets", "n", "_min", "_max")

    def __init__(self, alpha=0.01):
        if not 0.0 < alpha < 1.0:
            raise ValueError("alpha must be in (0, 1)")
        self.alpha = alpha
        self._gamma = (1.0 + alpha) / (1.0 - alpha)
        self._log_gamma = math.log(self._gamma)
        self._buckets = {}          # bucket index -> count
        self.n = 0
        self._min = math.inf
        self._max = -math.inf

    def _key(self, x):
        return int(math.ceil(math.log(x) / self._log_gamma))

    def add(self, x, weight=1):
        """Fold ``weight`` occurrences of a strictly positive value ``x``."""
        x = float(x)
        if x <= 0.0:
            raise ValueError("DDSketch handles strictly positive values only")
        if weight < 1:
            raise ValueError("weight must be a positive integer")
        k = self._key(x)
        self._buckets[k] = self._buckets.get(k, 0) + weight
        self.n += weight
        if x < self._min:
            self._min = x
        if x > self._max:
            self._max = x
        return self

    def quantile(self, q):
        """Estimate the ``q``-quantile (``0 <= q <= 1``) within relative error ``alpha``."""
        if self.n == 0:
            raise ValueError("cannot take a quantile of an empty sketch")
        if not 0.0 <= q <= 1.0:
            raise ValueError("q must be in [0, 1]")
        # rank of the desired element, 0-indexed, matching the sorted-array convention
        rank = q * (self.n - 1)
        target = int(math.floor(rank)) + 1     # count of elements up to and including it
        cumulative = 0
        for k in sorted(self._buckets):
            cumulative += self._buckets[k]
            if cumulative >= target:
                # value at the low edge of bucket k, scaled to the bucket midpoint region:
                # x_hat = 2 * gamma^k / (gamma + 1) keeps relative error <= alpha
                x_hat = 2.0 * (self._gamma ** k) / (self._gamma + 1.0)
                # clamp to the observed range so extreme quantiles never exceed it
                if x_hat < self._min:
                    return self._min
                if x_hat > self._max:
                    return self._max
                return x_hat
        return self._max

    def min(self):
        """Smallest value seen (exact)."""
        if self.n == 0:
            raise ValueError("empty sketch")
        return self._min

    def max(self):
        """Largest value seen (exact)."""
        if self.n == 0:
            raise ValueError("empty sketch")
        return self._max

    def merge(self, other):
        """Fold another sketch (same ``alpha``) into this one in place."""
        if not isinstance(other, DDSketch):
            raise TypeError("can only merge with another DDSketch")
        if abs(self.alpha - other.alpha) > 1e-15:
            raise ValueError("cannot merge sketches with different alpha")
        for k, c in other._buckets.items():
            self._buckets[k] = self._buckets.get(k, 0) + c
        self.n += other.n
        if other.n:
            self._min = min(self._min, other._min)
            self._max = max(self._max, other._max)
        return self

    def __add__(self, other):
        merged = DDSketch(self.alpha)
        merged.merge(self)
        merged.merge(other)
        return merged

    def num_buckets(self):
        """Number of populated buckets (the sketch's live memory footprint)."""
        return len(self._buckets)
