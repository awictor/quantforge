"""Extreme value theory: Hill tail index and peaks-over-threshold VaR / ES.

Tail risk from the extremes of a loss sample. The Hill estimator gauges the tail
index of a heavy (power-law) tail; the peaks-over-threshold approach fits a
Generalized Pareto distribution to exceedances over a high threshold and gives
closed-form tail VaR and expected shortfall. Losses are positive numbers. Pure
standard library.
"""

import math


def hill_estimator(losses, k):
    """Hill estimator of the tail index from the top ``k`` order statistics.

    ``xi = (1/k) sum_{i=1}^{k} ln(X_(n-i+1) / X_(n-k))`` -- the mean log-excess of
    the ``k`` largest losses over the ``(k+1)``-th. Estimates the shape ``xi`` of a
    heavy power-law tail (tail exponent ``alpha = 1/xi``); larger ``xi`` means a
    heavier tail.
    """
    xs = sorted(x for x in losses if x > 0)
    n = len(xs)
    if k < 1 or k >= n:
        raise ValueError("k must be in [1, n-1] with positive losses")
    top = xs[n - k:]
    threshold = xs[n - k - 1]
    if threshold <= 0:
        raise ValueError("threshold order statistic must be positive")
    return sum(math.log(x / threshold) for x in top) / k


def gpd_fit_pot(losses, threshold):
    """Fit a Generalized Pareto to peaks over ``threshold`` (method of moments).

    Returns ``(xi, beta, n_exceed, n_total)``: the GPD shape ``xi`` and scale
    ``beta`` matching the mean and variance of the exceedances ``X - threshold``,
    the number of exceedances, and the sample size. ``xi > 0`` is a heavy
    (Pareto-type) tail.
    """
    n = len(losses)
    if n == 0:
        raise ValueError("losses must be non-empty")
    exceed = [x - threshold for x in losses if x > threshold]
    ne = len(exceed)
    if ne < 2:
        raise ValueError("need at least two exceedances over the threshold")
    mean = sum(exceed) / ne
    var = sum((e - mean) ** 2 for e in exceed) / ne
    if var <= 0.0:
        raise ValueError("zero-variance exceedances")
    # Method-of-moments GPD: xi = 0.5 (1 - mean^2/var), beta = 0.5 mean (mean^2/var + 1).
    ratio = mean * mean / var
    xi = 0.5 * (1.0 - ratio)
    beta = 0.5 * mean * (ratio + 1.0)
    return xi, beta, ne, n


def gev_cdf(x, loc, scale, shape):
    """Generalized extreme value CDF (block-maxima limit distribution).

    ``F(x) = exp(-(1 + shape (x-loc)/scale)^{-1/shape})`` for ``shape != 0`` (Frechet
    ``shape > 0`` / Weibull ``shape < 0``), and the Gumbel limit
    ``exp(-e^{-(x-loc)/scale})`` as ``shape -> 0``. Defined where
    ``1 + shape (x-loc)/scale > 0``.
    """
    if scale <= 0:
        raise ValueError("scale must be positive")
    z = (x - loc) / scale
    if abs(shape) < 1e-8:
        return math.exp(-math.exp(-z))
    arg = 1.0 + shape * z
    if arg <= 0.0:
        return 0.0 if shape > 0 else 1.0
    return math.exp(-arg ** (-1.0 / shape))


def gev_return_level(period, loc, scale, shape):
    """Return level: the block maximum exceeded once per ``period`` blocks.

    Inverts :func:`gev_cdf` at ``p = 1 - 1/period``:

        level = loc + (scale/shape) [ (-ln(1 - 1/period))^{-shape} - 1 ],

    with the Gumbel limit ``loc - scale ln(-ln(1 - 1/period))`` as ``shape -> 0``.
    The T-block return level rises with the return period.
    """
    if period <= 1:
        raise ValueError("period must exceed 1")
    if scale <= 0:
        raise ValueError("scale must be positive")
    y = -math.log(1.0 - 1.0 / period)
    if abs(shape) < 1e-8:
        return loc - scale * math.log(y)
    return loc + (scale / shape) * (y ** (-shape) - 1.0)


def gev_fit_block_maxima(block_maxima):
    """Method-of-moments GEV fit assuming the Gumbel (shape = 0) limit.

    Fits location and scale of a Gumbel to the block maxima by moment matching:
    ``scale = std * sqrt(6)/pi``, ``loc = mean - gamma * scale`` (``gamma`` the
    Euler-Mascheroni constant). Returns ``(loc, scale, shape=0.0)`` -- a simple,
    robust baseline when the shape is not separately estimated.
    """
    n = len(block_maxima)
    if n < 2:
        raise ValueError("need at least two block maxima")
    mean = sum(block_maxima) / n
    var = sum((x - mean) ** 2 for x in block_maxima) / n
    if var <= 0.0:
        raise ValueError("zero-variance block maxima")
    scale = math.sqrt(var) * math.sqrt(6.0) / math.pi
    gamma = 0.5772156649015329
    loc = mean - gamma * scale
    return loc, scale, 0.0


def gpd_var(losses, threshold, confidence=0.99):
    """Peaks-over-threshold VaR from a fitted Generalized Pareto tail.

    ``VaR = threshold + (beta/xi) [ (n/N_u (1 - confidence))^{-xi} - 1 ]`` with
    ``N_u`` exceedances of ``n`` observations. Reduces to the exponential-tail
    limit as ``xi -> 0``. A positive loss quantile deep in the tail.
    """
    xi, beta, nu, n = gpd_fit_pot(losses, threshold)
    frac = n / nu * (1.0 - confidence)
    if abs(xi) < 1e-8:
        return threshold - beta * math.log(frac)
    return threshold + (beta / xi) * (frac ** (-xi) - 1.0)


def gpd_expected_shortfall(losses, threshold, confidence=0.99):
    """Peaks-over-threshold expected shortfall from the GPD tail.

    ``ES = VaR/(1 - xi) + (beta - xi threshold)/(1 - xi)`` for ``xi < 1``. Always
    at least the :func:`gpd_var`; finite only for a tail with ``xi < 1``.
    """
    xi, beta, _nu, _n = gpd_fit_pot(losses, threshold)
    if xi >= 1.0:
        raise ValueError("expected shortfall is infinite for xi >= 1")
    var = gpd_var(losses, threshold, confidence)
    return var / (1.0 - xi) + (beta - xi * threshold) / (1.0 - xi)
