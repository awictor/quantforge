"""Exponentially-weighted moving mean and variance (streaming).

An EWMA weights recent observations more than old ones, so it tracks a drifting mean or
a changing volatility without a fixed window. This is the online, single-series
accumulator (the RiskMetrics volatility recursion): each update blends the new point
into the running estimates by a decay ``lambda``,

    mean_t = lambda mean_{t-1} + (1-lambda) x_t,
    var_t  = lambda var_{t-1}  + (1-lambda) (x_t - mean_{t-1})^2.

The effective memory is about ``1/(1-lambda)`` points, so ``lambda = 0.94`` (RiskMetrics
daily) keeps roughly a month. Complements the batch ``ewma_covariance``. Pure standard
library.
"""

import math


class EWMAStats:
    """Streaming exponentially-weighted mean, variance and volatility.

    ``lam`` is the decay in ``(0, 1)`` (larger = longer memory). Feed points with
    :meth:`update`; read :attr:`mean`, :meth:`variance`, :meth:`std`. The first point
    seeds the mean; variance builds from the second. Effective window is about
    ``1 / (1 - lam)`` observations.
    """

    def __init__(self, lam=0.94, values=None):
        if not (0.0 < lam < 1.0):
            raise ValueError("lam must be in (0, 1)")
        self.lam = lam
        self.n = 0
        self.mean = 0.0
        self._var = 0.0
        if values is not None:
            for v in values:
                self.update(v)

    def update(self, x):
        """Incorporate one observation ``x`` and return the current mean."""
        self.n += 1
        if self.n == 1:
            self.mean = float(x)
            self._var = 0.0
            return self.mean
        prev_mean = self.mean
        # Variance uses the deviation from the *previous* mean (RiskMetrics form).
        dev = x - prev_mean
        self._var = self.lam * self._var + (1.0 - self.lam) * dev * dev
        self.mean = self.lam * prev_mean + (1.0 - self.lam) * x
        return self.mean

    def variance(self):
        """Current exponentially-weighted variance."""
        return self._var

    def std(self):
        """Current exponentially-weighted standard deviation (volatility)."""
        return math.sqrt(self._var) if self._var > 0 else 0.0


def ewma(values, lam=0.94):
    """Exponentially-weighted moving mean of a whole series (list output).

    Returns the running EWMA at each step, seeded from the first value. A convenience
    wrapper over :class:`EWMAStats` for batch use.
    """
    if not values:
        raise ValueError("need at least one value")
    acc = EWMAStats(lam=lam)
    return [acc.update(v) for v in values]
