"""Streaming (online) mean, variance, skewness, and kurtosis.

The batch skewness/kurtosis routines need the whole sample in memory and two passes;
this accumulator updates all four moments in one pass as data streams in, using the
numerically-stable central-moment recurrences (Terriberry's extension of Welford).
Two accumulators can also be *merged*, so partial results from parallel chunks combine
exactly into the whole -- the associativity that makes it a parallel reduction. Pure
standard library.
"""

import math


class RunningMoments:
    """One-pass accumulator for mean, variance, skewness, and excess kurtosis.

    Feed values with :meth:`update` (or a whole iterable to the constructor). Read
    :attr:`mean`, :meth:`variance`, :meth:`skewness`, :meth:`kurtosis` at any time.
    ``+`` merges two accumulators into one covering both samples, exactly.
    """

    def __init__(self, values=None):
        self.n = 0
        self.mean = 0.0
        self.M2 = 0.0        # sum of squared deviations
        self.M3 = 0.0        # third central moment sum
        self.M4 = 0.0        # fourth central moment sum
        if values is not None:
            for v in values:
                self.update(v)

    def update(self, x):
        """Incorporate one observation ``x`` (Terriberry's stable recurrence)."""
        n1 = self.n
        self.n += 1
        delta = x - self.mean
        delta_n = delta / self.n
        delta_n2 = delta_n * delta_n
        term1 = delta * delta_n * n1
        self.mean += delta_n
        self.M4 += (term1 * delta_n2 * (self.n * self.n - 3 * self.n + 3)
                    + 6 * delta_n2 * self.M2 - 4 * delta_n * self.M3)
        self.M3 += term1 * delta_n * (self.n - 2) - 3 * delta_n * self.M2
        self.M2 += term1
        return self

    def variance(self, sample=True):
        """Variance; sample (``n-1``) by default, population if ``sample=False``."""
        if self.n < 2:
            return 0.0
        return self.M2 / (self.n - 1 if sample else self.n)

    def std(self, sample=True):
        """Standard deviation (square root of :meth:`variance`)."""
        return math.sqrt(self.variance(sample))

    def skewness(self):
        """Population skewness ``sqrt(n) M3 / M2^{1.5}`` (0 for a symmetric sample)."""
        if self.n < 2 or self.M2 == 0.0:
            return 0.0
        return math.sqrt(self.n) * self.M3 / (self.M2 ** 1.5)

    def kurtosis(self):
        """Excess kurtosis ``n M4 / M2^2 - 3`` (0 for a normal sample)."""
        if self.n < 2 or self.M2 == 0.0:
            return 0.0
        return self.n * self.M4 / (self.M2 * self.M2) - 3.0

    def __add__(self, other):
        """Merge two accumulators exactly (Chan et al. parallel combination)."""
        if not isinstance(other, RunningMoments):
            return NotImplemented
        a, b = self, other
        if a.n == 0:
            return _clone(b)
        if b.n == 0:
            return _clone(a)
        n = a.n + b.n
        delta = b.mean - a.mean
        d2 = delta * delta
        d3 = d2 * delta
        d4 = d2 * d2
        r = RunningMoments()
        r.n = n
        r.mean = a.mean + delta * b.n / n
        r.M2 = a.M2 + b.M2 + d2 * a.n * b.n / n
        r.M3 = (a.M3 + b.M3
                + d3 * a.n * b.n * (a.n - b.n) / (n * n)
                + 3 * delta * (a.n * b.M2 - b.n * a.M2) / n)
        r.M4 = (a.M4 + b.M4
                + d4 * a.n * b.n * (a.n * a.n - a.n * b.n + b.n * b.n) / (n ** 3)
                + 6 * d2 * (a.n * a.n * b.M2 + b.n * b.n * a.M2) / (n * n)
                + 4 * delta * (a.n * b.M3 - b.n * a.M3) / n)
        return r


def _clone(src):
    r = RunningMoments()
    r.n, r.mean, r.M2, r.M3, r.M4 = src.n, src.mean, src.M2, src.M3, src.M4
    return r
