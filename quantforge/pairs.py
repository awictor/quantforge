"""Pairs trading: hedge ratio, spread, and Ornstein-Uhlenbeck half-life.

For a mean-reverting pair, the hedge ratio is the OLS slope of one leg on the
other; the spread is the residual. Fitting an AR(1) to the spread gives the
mean-reversion speed and the half-life of a shock. This module computes the hedge
ratio, spread, and OU half-life, plus the spread z-score for entry/exit. Pure
standard library.
"""

import math


def pairs_hedge_ratio(y, x):
    """OLS hedge ratio (slope) of ``y`` on ``x`` through the mean.

    ``beta = cov(x, y) / var(x)`` -- the number of units of ``x`` to short against
    one unit of ``y`` so the spread ``y - beta x`` is mean-reverting. Recovers the
    true beta on a linear relationship.
    """
    n = len(y)
    if len(x) != n:
        raise ValueError("series must have equal length")
    if n < 2:
        raise ValueError("need at least two observations")
    mx = sum(x) / n
    my = sum(y) / n
    var_x = sum((xi - mx) ** 2 for xi in x)
    if var_x <= 0.0:
        raise ValueError("x must have positive variance")
    cov = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    return cov / var_x


def spread_series(y, x, beta=None):
    """Spread ``y - beta x`` (hedge-ratio residual).

    ``beta`` defaults to the :func:`hedge_ratio`. The series a pairs trade bets
    reverts to its mean.
    """
    if beta is None:
        beta = pairs_hedge_ratio(y, x)
    return [y[i] - beta * x[i] for i in range(len(y))]


def ou_half_life(spread):
    """Mean-reversion half-life from an AR(1) fit to the spread.

    Regresses ``delta_t = a + b * spread_{t-1}`` (the discretized OU); the
    mean-reversion speed is ``kappa = -b`` and the half-life is ``ln(2)/kappa``.
    Positive and finite only for a mean-reverting (``-1 < b < 0``) spread; raises
    otherwise.
    """
    n = len(spread)
    if n < 3:
        raise ValueError("need at least three spread observations")
    lag = spread[:-1]
    delta = [spread[i + 1] - spread[i] for i in range(n - 1)]
    m_lag = sum(lag) / len(lag)
    m_del = sum(delta) / len(delta)
    var_lag = sum((l - m_lag) ** 2 for l in lag)
    if var_lag <= 0.0:
        raise ValueError("spread lag has zero variance")
    b = sum((lag[i] - m_lag) * (delta[i] - m_del) for i in range(len(lag))) / var_lag
    kappa = -b
    if kappa <= 0.0:
        raise ValueError("spread is not mean-reverting (non-positive kappa)")
    return math.log(2.0) / kappa


def spread_zscore(spread, window=None):
    """Latest spread z-score against its mean and std (full history or a window).

    ``(spread[-1] - mean) / std`` over the last ``window`` points (all if ``None``).
    The pairs-trade entry signal: large magnitude means the spread is stretched.
    """
    s = spread if window is None else spread[-window:]
    n = len(s)
    if n < 2:
        raise ValueError("need at least two observations")
    m = sum(s) / n
    var = sum((x - m) ** 2 for x in s) / (n - 1)
    sd = var ** 0.5
    if sd == 0.0:
        return 0.0
    return (s[-1] - m) / sd
