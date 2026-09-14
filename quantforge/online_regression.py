"""Online (single-pass) simple linear regression.

Fits ``y = slope * x + intercept`` by ordinary least squares over a stream of ``(x, y)``
pairs in one pass and constant memory, using the numerically stable Welford co-moment
recursion (the same update as :class:`RunningCovariance`). Reads slope, intercept,
correlation, and R-squared at any point; predicts for new ``x``; and two accumulators
merge exactly, so chunks fitted in parallel combine into the whole-sample fit. Pure
standard library.
"""

import math


class RunningRegression:
    """Streaming OLS fit of ``y`` on ``x`` (Welford co-moment, constant memory).

    ``update(x, y)`` folds in one pair; ``slope()``, ``intercept()``, ``correlation()``,
    ``r_squared()``, and ``predict(x)`` read the current fit at any time. Needs at least
    two points with variation in ``x`` before a slope is defined. Two accumulators combine
    with ``+`` via the parallel-merge formula, giving the exact same fit as folding every
    pair into one.
    """

    __slots__ = ("n", "mean_x", "mean_y", "c_xy", "m2_x", "m2_y")

    def __init__(self, xs=None, ys=None):
        self.n = 0
        self.mean_x = 0.0
        self.mean_y = 0.0
        self.c_xy = 0.0          # sum of co-moments (x-mx)(y-my)
        self.m2_x = 0.0          # sum of squared deviations of x
        self.m2_y = 0.0          # sum of squared deviations of y
        if xs is not None:
            if ys is None or len(xs) != len(ys):
                raise ValueError("xs and ys must have equal length")
            for x, y in zip(xs, ys):
                self.update(x, y)

    def update(self, x, y):
        """Fold one ``(x, y)`` pair into the running fit."""
        x = float(x)
        y = float(y)
        self.n += 1
        n = self.n
        dx = x - self.mean_x
        dy = y - self.mean_y
        self.mean_x += dx / n
        self.mean_y += dy / n
        # co-moment and second moments use the post-update mean for one factor
        self.c_xy += dx * (y - self.mean_y)
        self.m2_x += dx * (x - self.mean_x)
        self.m2_y += dy * (y - self.mean_y)
        return self

    def slope(self):
        """OLS slope ``cov(x, y) / var(x)``; raises if fewer than 2 points or x has no spread."""
        if self.n < 2:
            raise ValueError("need at least two points for a slope")
        if self.m2_x <= 0.0:
            raise ValueError("x has zero variance; slope undefined")
        return self.c_xy / self.m2_x

    def intercept(self):
        """OLS intercept ``mean_y - slope * mean_x``."""
        return self.mean_y - self.slope() * self.mean_x

    def predict(self, x):
        """Predicted ``y`` for a new ``x`` from the current fit."""
        return self.intercept() + self.slope() * float(x)

    def correlation(self):
        """Pearson correlation of the pairs seen so far, in ``[-1, 1]``."""
        if self.n < 2:
            raise ValueError("need at least two points for a correlation")
        denom = self.m2_x * self.m2_y
        if denom <= 0.0:
            raise ValueError("a variable has zero variance; correlation undefined")
        r = self.c_xy / math.sqrt(denom)
        if r > 1.0:
            return 1.0
        if r < -1.0:
            return -1.0
        return r

    def r_squared(self):
        """Coefficient of determination; equals the squared correlation in simple OLS."""
        r = self.correlation()
        return r * r

    def __add__(self, other):
        """Merge two accumulators (parallel-merge), giving the combined-sample fit."""
        if not isinstance(other, RunningRegression):
            return NotImplemented
        if self.n == 0:
            return other._copy()
        if other.n == 0:
            return self._copy()
        merged = RunningRegression()
        na, nb = self.n, other.n
        n = na + nb
        delta_x = other.mean_x - self.mean_x
        delta_y = other.mean_y - self.mean_y
        merged.n = n
        merged.mean_x = self.mean_x + delta_x * nb / n
        merged.mean_y = self.mean_y + delta_y * nb / n
        f = na * nb / n
        merged.c_xy = self.c_xy + other.c_xy + delta_x * delta_y * f
        merged.m2_x = self.m2_x + other.m2_x + delta_x * delta_x * f
        merged.m2_y = self.m2_y + other.m2_y + delta_y * delta_y * f
        return merged

    def _copy(self):
        c = RunningRegression()
        c.n = self.n
        c.mean_x = self.mean_x
        c.mean_y = self.mean_y
        c.c_xy = self.c_xy
        c.m2_x = self.m2_x
        c.m2_y = self.m2_y
        return c
