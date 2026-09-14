"""Online (single-pass) covariance and correlation.

Tracks the covariance and Pearson correlation of a stream of ``(x, y)`` pairs in one pass
and constant memory, using the numerically stable Welford co-moment update. Two
accumulators merge exactly (parallel/chunked data), like :class:`RunningMoments`. Pure
standard library.
"""

import math


class RunningCovariance:
    """Streaming covariance/correlation of paired observations (Welford co-moment).

    ``update(x, y)`` folds in one pair; ``covariance()`` and ``correlation()`` read the
    current estimate at any time. Sample (``ddof=1``) covariance by default. Two
    accumulators combine with ``+`` using the parallel-merge formula.
    """

    __slots__ = ("n", "mean_x", "mean_y", "c_xy", "m2_x", "m2_y")

    def __init__(self, xs=None, ys=None):
        self.n = 0
        self.mean_x = 0.0
        self.mean_y = 0.0
        self.c_xy = 0.0          # sum of co-moments
        self.m2_x = 0.0          # sum of squared deviations of x
        self.m2_y = 0.0
        if xs is not None:
            if ys is None or len(xs) != len(ys):
                raise ValueError("xs and ys must be equal-length")
            for x, y in zip(xs, ys):
                self.update(x, y)

    def update(self, x, y):
        """Fold one ``(x, y)`` pair into the running estimates."""
        self.n += 1
        dx = x - self.mean_x
        self.mean_x += dx / self.n
        dy = y - self.mean_y
        self.mean_y += dy / self.n
        self.c_xy += dx * (y - self.mean_y)
        self.m2_x += dx * (x - self.mean_x)
        self.m2_y += dy * (y - self.mean_y)
        return self

    def covariance(self, ddof=1):
        """Covariance of the pairs seen so far (``ddof=1`` sample, ``0`` population)."""
        if self.n - ddof <= 0:
            raise ValueError("not enough data for the requested ddof")
        return self.c_xy / (self.n - ddof)

    def variance_x(self, ddof=1):
        """Variance of the ``x`` marginal."""
        if self.n - ddof <= 0:
            raise ValueError("not enough data for the requested ddof")
        return self.m2_x / (self.n - ddof)

    def variance_y(self, ddof=1):
        """Variance of the ``y`` marginal."""
        if self.n - ddof <= 0:
            raise ValueError("not enough data for the requested ddof")
        return self.m2_y / (self.n - ddof)

    def correlation(self):
        """Pearson correlation of the pairs so far, in ``[-1, 1]``.

        Raises if either marginal has zero variance.
        """
        if self.m2_x <= 0.0 or self.m2_y <= 0.0:
            raise ValueError("a marginal has zero variance")
        return self.c_xy / math.sqrt(self.m2_x * self.m2_y)

    @property
    def means(self):
        """Current ``(mean_x, mean_y)``."""
        return self.mean_x, self.mean_y

    def __add__(self, other):
        """Merge two accumulators exactly (Chan parallel co-moment formula)."""
        if self.n == 0:
            return _copy(other)
        if other.n == 0:
            return _copy(self)
        merged = RunningCovariance()
        n = self.n + other.n
        dx = other.mean_x - self.mean_x
        dy = other.mean_y - self.mean_y
        merged.n = n
        merged.mean_x = self.mean_x + dx * other.n / n
        merged.mean_y = self.mean_y + dy * other.n / n
        f = self.n * other.n / n
        merged.c_xy = self.c_xy + other.c_xy + dx * dy * f
        merged.m2_x = self.m2_x + other.m2_x + dx * dx * f
        merged.m2_y = self.m2_y + other.m2_y + dy * dy * f
        return merged


def _copy(acc):
    out = RunningCovariance()
    out.n = acc.n
    out.mean_x = acc.mean_x
    out.mean_y = acc.mean_y
    out.c_xy = acc.c_xy
    out.m2_x = acc.m2_x
    out.m2_y = acc.m2_y
    return out
