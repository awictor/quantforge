"""Trend and momentum signals: moving averages, MACD, RSI, z-score.

Classic technical/systematic signals computed from a price or return series:
simple and exponential moving averages, the MACD oscillator, the relative
strength index, a rolling z-score, and the time-series-momentum sign. Pure
standard library.
"""


def sma(series, window):
    """Simple moving average over each trailing ``window`` (list, one per position).

    Returns ``len(series) - window + 1`` values, each the mean of that window.
    """
    n = len(series)
    if window < 1 or window > n:
        raise ValueError("window must be in [1, len(series)]")
    return [sum(series[i:i + window]) / window for i in range(n - window + 1)]


def ema(series, span):
    """Exponential moving average with smoothing ``alpha = 2/(span+1)``.

    Recursive ``e_t = alpha x_t + (1 - alpha) e_{t-1}`` seeded at the first point.
    Reacts faster than the :func:`sma` of the same length (less lag).
    """
    n = len(series)
    if span < 1 or n == 0:
        raise ValueError("span >= 1 and non-empty series required")
    alpha = 2.0 / (span + 1.0)
    out = [series[0]]
    for x in series[1:]:
        out.append(alpha * x + (1.0 - alpha) * out[-1])
    return out


def macd(series, fast=12, slow=26, signal=9):
    """MACD line, signal line, and histogram.

    MACD line = ``EMA(fast) - EMA(slow)``; signal line = ``EMA(signal)`` of the
    MACD line; histogram = MACD - signal. Returns ``(macd_line, signal_line,
    histogram)`` aligned to the series. Positive MACD indicates the fast average
    above the slow (up-momentum).
    """
    if fast >= slow:
        raise ValueError("fast span must be less than slow span")
    ef, es = ema(series, fast), ema(series, slow)
    macd_line = [ef[i] - es[i] for i in range(len(series))]
    signal_line = ema(macd_line, signal)
    hist = [macd_line[i] - signal_line[i] for i in range(len(series))]
    return macd_line, signal_line, hist


def rsi(series, window=14):
    """Relative strength index over a trailing ``window`` (Wilder's smoothing).

    ``RSI = 100 - 100/(1 + avg_gain/avg_loss)`` in ``[0, 100]``. Above 70 is
    conventionally overbought, below 30 oversold; near 100 in a strong uptrend.
    Returns one value per position from index ``window`` on.
    """
    n = len(series)
    if window < 1 or window >= n:
        raise ValueError("window must be in [1, len(series) - 1]")
    gains = [max(series[i] - series[i - 1], 0.0) for i in range(1, n)]
    losses = [max(series[i - 1] - series[i], 0.0) for i in range(1, n)]
    avg_gain = sum(gains[:window]) / window
    avg_loss = sum(losses[:window]) / window
    out = []
    for i in range(window, n):
        if i > window:
            avg_gain = (avg_gain * (window - 1) + gains[i - 1]) / window
            avg_loss = (avg_loss * (window - 1) + losses[i - 1]) / window
        if avg_loss == 0.0:
            out.append(100.0)
        else:
            rs = avg_gain / avg_loss
            out.append(100.0 - 100.0 / (1.0 + rs))
    return out


def rolling_zscore(series, window):
    """Rolling z-score ``(x_t - mean) / std`` over each trailing ``window``.

    Standardizes the latest point against its window; a mean-reversion / breakout
    signal. Windows with zero variance yield 0. One value per position from index
    ``window - 1`` on.
    """
    n = len(series)
    if window < 2 or window > n:
        raise ValueError("window must be in [2, len(series)]")
    out = []
    for end in range(window, n + 1):
        w = series[end - window:end]
        m = sum(w) / window
        var = sum((x - m) ** 2 for x in w) / (window - 1)
        sd = var ** 0.5
        out.append(0.0 if sd == 0.0 else (w[-1] - m) / sd)
    return out


def bollinger_bands(series, window=20, num_std=2.0):
    """Bollinger bands: ``(lower, middle, upper)`` lists over a trailing window.

    Middle is the :func:`sma`; the bands are ``middle +/- num_std * rolling std``.
    Price closing above the upper / below the lower band flags stretched moves.
    Returns three aligned lists, one value per window position.
    """
    n = len(series)
    if window < 2 or window > n:
        raise ValueError("window must be in [2, len(series)]")
    lower, middle, upper = [], [], []
    for end in range(window, n + 1):
        w = series[end - window:end]
        m = sum(w) / window
        var = sum((x - m) ** 2 for x in w) / window
        sd = var ** 0.5
        middle.append(m)
        lower.append(m - num_std * sd)
        upper.append(m + num_std * sd)
    return lower, middle, upper


def average_true_range(highs, lows, closes, window=14):
    """Average true range (Wilder): mean of the true range over a trailing window.

    True range at ``t`` is ``max(high-low, |high-prev_close|, |low-prev_close|)``;
    ATR smooths it with Wilder's moving average. A non-negative volatility measure
    in price units. Returns one value per position from index ``window`` on.
    """
    n = len(closes)
    if not (len(highs) == len(lows) == n):
        raise ValueError("highs, lows, closes must have equal length")
    if window < 1 or window >= n:
        raise ValueError("window must be in [1, len(closes) - 1]")
    tr = []
    for t in range(1, n):
        tr.append(max(highs[t] - lows[t], abs(highs[t] - closes[t - 1]),
                      abs(lows[t] - closes[t - 1])))
    atr = sum(tr[:window]) / window
    out = [atr]
    for t in range(window, len(tr)):
        atr = (atr * (window - 1) + tr[t]) / window
        out.append(atr)
    return out


def donchian_channel(highs, lows, window=20):
    """Donchian channel: rolling ``(lowest_low, highest_high)`` over a window.

    The channel a breakout system trades: a close above the prior highest high is
    a long breakout, below the lowest low a short. Returns ``(lower, upper)`` lists
    with ``upper >= lower`` at every position.
    """
    n = len(highs)
    if len(lows) != n:
        raise ValueError("highs and lows must have equal length")
    if window < 1 or window > n:
        raise ValueError("window must be in [1, len(highs)]")
    lower, upper = [], []
    for end in range(window, n + 1):
        lower.append(min(lows[end - window:end]))
        upper.append(max(highs[end - window:end]))
    return lower, upper


def time_series_momentum(prices, lookback):
    """Sign of the trailing ``lookback``-period return: +1 up, -1 down, 0 flat.

    The time-series-momentum signal (Moskowitz-Ooi-Pedersen): go long after a
    positive past return, short after a negative one. Returns one signal per
    position from index ``lookback`` on.
    """
    n = len(prices)
    if lookback < 1 or lookback >= n:
        raise ValueError("lookback must be in [1, len(prices) - 1]")
    out = []
    for t in range(lookback, n):
        past = prices[t] - prices[t - lookback]
        out.append(1 if past > 0 else (-1 if past < 0 else 0))
    return out
