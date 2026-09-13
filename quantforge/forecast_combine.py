"""Combining forecasts: simple, inverse-MSE, and variance-minimizing weights.

Averaging competing forecasts usually beats any single one. Three schemes for the
combination weights, given each model's forecast-error series:

- ``simple_average_forecast`` -- equal weights,
- ``inverse_mse_weights`` -- weights proportional to ``1 / MSE``, ignoring the
  cross-model error correlation,
- ``optimal_combination_weights`` -- Bates-Granger minimum-variance weights
  ``w = Sigma^{-1} 1 / (1' Sigma^{-1} 1)`` from the error covariance ``Sigma``,
  which is optimal when the errors are correlated.

Weights sum to one in every scheme. Pure standard library.
"""

from .portopt import _invert


def _mse(errors):
    n = len(errors)
    return sum(e * e for e in errors) / n


def simple_average_forecast(forecasts):
    """Equal-weight combination of aligned forecast series.

    ``forecasts`` is a list of series (one per model); returns the point-wise mean.
    """
    m = len(forecasts)
    if m == 0:
        raise ValueError("need at least one forecast series")
    n = len(forecasts[0])
    if any(len(f) != n for f in forecasts):
        raise ValueError("forecast series must be equal length")
    return [sum(forecasts[j][i] for j in range(m)) / m for i in range(n)]


def inverse_mse_weights(error_series):
    """Combination weights proportional to ``1 / MSE`` of each model's errors.

    ``error_series`` is a list of forecast-error series. A more accurate model (lower
    MSE) gets more weight; weights sum to one. Ignores cross-model error correlation.
    """
    m = len(error_series)
    if m == 0:
        raise ValueError("need at least one error series")
    inv = [1.0 / _mse(e) if _mse(e) > 0 else float("inf") for e in error_series]
    total = sum(inv)
    if total == float("inf"):
        # A perfect model dominates; put all weight there.
        return [1.0 if x == float("inf") else 0.0 for x in inv]
    return [x / total for x in inv]


def optimal_combination_weights(error_series):
    """Bates-Granger minimum-variance combination weights.

    ``w = Sigma^{-1} 1 / (1' Sigma^{-1} 1)`` where ``Sigma`` is the covariance of the
    models' forecast errors. Minimizes the variance of the combined error and can put
    negative weight on a model that hedges another's errors. Weights sum to one.
    """
    m = len(error_series)
    if m == 0:
        raise ValueError("need at least one error series")
    n = len(error_series[0])
    if n < 2 or any(len(e) != n for e in error_series):
        raise ValueError("error series must be equal length, >= 2 points")
    means = [sum(e) / n for e in error_series]
    cov = [[sum((error_series[i][t] - means[i]) * (error_series[j][t] - means[j])
                for t in range(n)) / n for j in range(m)] for i in range(m)]
    # Add the mean-outer-product so Sigma is the raw second moment (MSE on diagonal),
    # which is what Bates-Granger minimizes for biased forecasts too.
    for i in range(m):
        for j in range(m):
            cov[i][j] += means[i] * means[j]
    inv = _invert(cov)
    ones = [1.0] * m
    inv_ones = [sum(inv[i][j] * ones[j] for j in range(m)) for i in range(m)]
    denom = sum(inv_ones)
    if denom == 0.0:
        raise ValueError("degenerate error covariance")
    return [w / denom for w in inv_ones]


def combine(forecasts, weights):
    """Weighted combination of aligned forecast series by ``weights``."""
    m = len(forecasts)
    if m == 0 or len(weights) != m:
        raise ValueError("weights must align with the forecast series")
    n = len(forecasts[0])
    if any(len(f) != n for f in forecasts):
        raise ValueError("forecast series must be equal length")
    return [sum(weights[j] * forecasts[j][i] for j in range(m)) for i in range(n)]
