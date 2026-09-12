"""Exponential smoothing: Holt (trend) and Holt-Winters (seasonal).

Exponential smoothing forecasts a series as a weighted blend of its own recent
level, trend, and seasonal components, with geometrically decaying weights.

  * ``holt_linear`` -- double exponential smoothing with a level ``l`` and trend
    ``b``:  l_t = a y_t + (1-a)(l_{t-1}+b_{t-1});  b_t = g(l_t-l_{t-1})+(1-g)b_{t-1}.
    The h-step forecast is ``l_t + h b_t`` -- a straight line.
  * ``holt_winters_add`` -- adds an additive seasonal component of period ``m``:
    the forecast is ``l_t + h b_t + s_{t+h-m}``, capturing a repeating pattern.

Both are causal, O(n) recursions. Pure standard library.
"""


def holt_linear(y, alpha, beta, horizon=1):
    """Holt's linear-trend (double) exponential smoothing.

    Parameters
    ----------
    y : sequence of float
        The series (length >= 2).
    alpha, beta : float
        Level and trend smoothing parameters in [0, 1].
    horizon : int
        Forecast horizon (number of steps beyond the last observation).

    Returns
    -------
    (level, trend, forecast) : (float, float, list[float])
        Final level and trend, and the ``horizon``-step-ahead forecasts
        ``l + h*b`` (a straight line in ``h``).
    """
    n = len(y)
    if n < 2:
        raise ValueError("need at least 2 observations")
    if not (0.0 <= alpha <= 1.0 and 0.0 <= beta <= 1.0):
        raise ValueError("alpha and beta must be in [0, 1]")
    if horizon < 1:
        raise ValueError("horizon must be >= 1")

    level = y[0]
    trend = y[1] - y[0]
    for t in range(1, n):
        prev_level = level
        level = alpha * y[t] + (1.0 - alpha) * (level + trend)
        trend = beta * (level - prev_level) + (1.0 - beta) * trend
    forecast = [level + (h + 1) * trend for h in range(horizon)]
    return level, trend, forecast


def holt_winters_add(y, alpha, beta, gamma, period, horizon=1):
    """Additive Holt-Winters (triple) exponential smoothing.

    Parameters
    ----------
    y : sequence of float
        The series; length must exceed ``2 * period``.
    alpha, beta, gamma : float
        Level, trend, and seasonal smoothing parameters in [0, 1].
    period : int
        Season length ``m`` (>= 2).
    horizon : int
        Forecast horizon.

    Returns
    -------
    (level, trend, seasonals, forecast) : (float, float, list[float], list[float])
        Final level, trend, the last ``period`` seasonal factors, and the
        ``horizon``-step forecasts ``l + h*b + s[(h-1) mod m]``.
    """
    n = len(y)
    m = period
    if m < 2:
        raise ValueError("period must be >= 2")
    if n < 2 * m:
        raise ValueError("need more than two full seasons of data")
    if not all(0.0 <= p <= 1.0 for p in (alpha, beta, gamma)):
        raise ValueError("alpha, beta, gamma must be in [0, 1]")
    if horizon < 1:
        raise ValueError("horizon must be >= 1")

    # Initialize: level = mean of first season; trend = avg season-over-season
    # slope; seasonals = first-season deviations from that level.
    first = sum(y[:m]) / m
    second = sum(y[m:2 * m]) / m
    level = first
    trend = (second - first) / m
    seasonals = [y[i] - first for i in range(m)]

    for t in range(m, n):
        prev_level = level
        s_idx = (t - m) % m
        level = alpha * (y[t] - seasonals[s_idx]) + (1.0 - alpha) * (level + trend)
        trend = beta * (level - prev_level) + (1.0 - beta) * trend
        seasonals[t % m] = (gamma * (y[t] - level)
                            + (1.0 - gamma) * seasonals[s_idx])

    forecast = []
    for h in range(1, horizon + 1):
        s = seasonals[(n + h - 1) % m]
        forecast.append(level + h * trend + s)
    return level, trend, list(seasonals), forecast
