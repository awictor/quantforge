"""Heteroskedasticity- and autocorrelation-consistent (HAC) variance.

The sample variance understates the variance of a mean (or the long-run variance
of a series) when observations are autocorrelated. The Newey-West (1987) estimator
corrects for this by adding Bartlett-weighted autocovariances up to a lag ``L``:

    sigma_LR^2 = gamma_0 + 2 * sum_{k=1}^{L} (1 - k/(L+1)) * gamma_k,

where ``gamma_k`` is the lag-``k`` autocovariance. The Bartlett (triangular)
weights ``1 - k/(L+1)`` are exactly what guarantees the estimate is
non-negative -- an unweighted truncated sum can go negative. This is the standard
long-run variance behind HAC standard errors and spectral-density-at-frequency-
zero estimates. Pure standard library.
"""

import math


def autocovariance(x, lag):
    """Biased (divisor ``n``) sample autocovariance of ``x`` at ``lag``.

    ``gamma_k = (1/n) sum_{t=k}^{n-1} (x_t - xbar)(x_{t-k} - xbar)``. The divisor
    ``n`` (not ``n - k``) is what makes the autocovariance sequence positive
    semidefinite, which the Newey-West weighting relies on.
    """
    n = len(x)
    if lag < 0:
        raise ValueError("lag must be non-negative")
    if lag >= n:
        raise ValueError("lag must be smaller than the series length")
    mean = sum(x) / n
    acc = 0.0
    for t in range(lag, n):
        acc += (x[t] - mean) * (x[t - lag] - mean)
    return acc / n


def autocorrelation(x, lag):
    """Sample autocorrelation ``rho_k = gamma_k / gamma_0`` at ``lag``."""
    g0 = autocovariance(x, 0)
    if g0 == 0.0:
        raise ValueError("series has zero variance")
    return autocovariance(x, lag) / g0


def newey_west_variance(x, lags):
    """Newey-West long-run variance of a series (Bartlett-weighted HAC).

    Parameters
    ----------
    x : sequence of float
        The series (e.g. demeaned returns or a moment condition).
    lags : int
        Truncation lag ``L`` (>= 0). ``lags = 0`` reduces to the sample variance.

    Returns
    -------
    float
        The long-run variance estimate. Always non-negative thanks to the Bartlett
        weights.
    """
    n = len(x)
    if lags < 0:
        raise ValueError("lags must be non-negative")
    if lags >= n:
        raise ValueError("lags must be smaller than the series length")
    total = autocovariance(x, 0)
    for k in range(1, lags + 1):
        weight = 1.0 - k / (lags + 1.0)
        total += 2.0 * weight * autocovariance(x, k)
    return max(total, 0.0)


def newey_west_mean_se(x, lags):
    """HAC (Newey-West) standard error of the sample mean of ``x``.

    ``se = sqrt(sigma_LR^2 / n)`` with ``sigma_LR^2`` the Newey-West long-run
    variance. For serially uncorrelated data this matches the usual
    ``sqrt(var / n)``; positive autocorrelation inflates it.
    """
    n = len(x)
    return math.sqrt(newey_west_variance(x, lags) / n)
