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

Pure standard library.
"""


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
