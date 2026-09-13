"""Forecast-accuracy metrics.

Standard scale-dependent and scale-free errors for comparing forecasts against
realized values:

  * ``mae``   -- mean absolute error (same units as the data).
  * ``rmse``  -- root mean squared error (penalizes large misses more).
  * ``mape``  -- mean absolute percentage error (scale-free, undefined at zero
    actuals).
  * ``smape`` -- symmetric MAPE, bounded in [0, 2] (200%), robust to small actuals.
  * ``mase``  -- mean absolute scaled error: MAE divided by the in-sample MAE of a
    naive (random-walk) forecast. ``< 1`` beats naive, ``= 1`` matches it.
  * ``theil_u2`` -- Theil's U2: forecast RMSE over the no-change naive RMSE.
    ``< 1`` beats naive, ``= 1`` matches, ``> 1`` worse.
  * ``theil_u1`` -- Theil's U1 inequality coefficient, bounded in ``[0, 1]``.

Pure standard library.
"""

import math


def _check(actual, forecast):
    if len(actual) != len(forecast):
        raise ValueError("actual and forecast must have the same length")
    if len(actual) == 0:
        raise ValueError("need at least one point")


def mae(actual, forecast):
    """Mean absolute error."""
    _check(actual, forecast)
    n = len(actual)
    return sum(abs(actual[i] - forecast[i]) for i in range(n)) / n


def rmse(actual, forecast):
    """Root mean squared error."""
    _check(actual, forecast)
    n = len(actual)
    return (sum((actual[i] - forecast[i]) ** 2 for i in range(n)) / n) ** 0.5


def mape(actual, forecast):
    """Mean absolute percentage error (as a fraction; 0.1 = 10%).

    Raises if any actual value is zero (the percentage is undefined there).
    """
    _check(actual, forecast)
    n = len(actual)
    if any(a == 0.0 for a in actual):
        raise ValueError("mape undefined when an actual value is zero")
    return sum(abs((actual[i] - forecast[i]) / actual[i]) for i in range(n)) / n


def smape(actual, forecast):
    """Symmetric mean absolute percentage error, in [0, 2].

    ``mean( |a - f| / ((|a| + |f|) / 2) )``. Terms with both ``a`` and ``f`` zero
    contribute 0. Robust to small actuals and bounded, unlike plain MAPE.
    """
    _check(actual, forecast)
    n = len(actual)
    total = 0.0
    for i in range(n):
        denom = (abs(actual[i]) + abs(forecast[i])) / 2.0
        if denom > 0.0:
            total += abs(actual[i] - forecast[i]) / denom
    return total / n


def mase(actual, forecast, train=None, season=1):
    """Mean absolute scaled error.

    Scales the forecast MAE by the in-sample MAE of a seasonal-naive forecast
    (differences at lag ``season``) computed on ``train`` (defaults to ``actual``).
    ``< 1`` means the forecast beats naive; ``= 1`` matches it. Raises if the naive
    benchmark has zero error (a perfectly predictable training series).
    """
    _check(actual, forecast)
    hist = list(actual if train is None else train)
    if len(hist) <= season:
        raise ValueError("train shorter than the seasonal lag")
    naive_mae = sum(abs(hist[i] - hist[i - season])
                    for i in range(season, len(hist))) / (len(hist) - season)
    if naive_mae == 0.0:
        raise ValueError("naive benchmark has zero error; MASE undefined")
    return mae(actual, forecast) / naive_mae


def theil_u2(actual, forecast, last_actual=None):
    """Theil's U2 statistic: forecast RMSE over the no-change naive RMSE.

    The naive forecast for period ``t`` is the previous actual ``actual[t-1]`` (or
    ``last_actual`` for the first point, if given). ``U2 < 1`` means the forecast
    beats the random walk, ``= 1`` matches it, ``> 1`` is worse. Needs at least two
    points (or one point plus ``last_actual``).
    """
    _check(actual, forecast)
    n = len(actual)
    start = 1 if last_actual is None else 0
    prev0 = actual[0] if last_actual is None else last_actual
    num = 0.0
    den = 0.0
    prev = prev0
    for t in range(start, n):
        num += (forecast[t] - actual[t]) ** 2
        den += (prev - actual[t]) ** 2
        prev = actual[t]
    if den <= 0.0:
        raise ValueError("naive (no-change) forecast has zero error; U2 undefined")
    return math.sqrt(num / den)


def theil_u1(actual, forecast):
    """Theil's U1 inequality coefficient, bounded in ``[0, 1]``.

    ``U1 = RMSE / (rms(actual) + rms(forecast))``. Zero for a perfect forecast and
    one in the worst case; scale-free and symmetric. Distinct from
    :func:`theil_u2`, which benchmarks against the naive forecast.
    """
    _check(actual, forecast)
    n = len(actual)
    mse = sum((forecast[t] - actual[t]) ** 2 for t in range(n)) / n
    rms_a = math.sqrt(sum(a * a for a in actual) / n)
    rms_f = math.sqrt(sum(f * f for f in forecast) / n)
    denom = rms_a + rms_f
    if denom <= 0.0:
        return 0.0
    return math.sqrt(mse) / denom
