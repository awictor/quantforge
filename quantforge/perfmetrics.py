"""Performance and drawdown statistics for a return series.

Standard track-record measures computed from a sequence of periodic returns
(simple, not log): the Sharpe and Sortino ratios (annualized), maximum
drawdown of the cumulative-return curve, the Calmar ratio, and hit rate /
profit factor. All pure standard library.
"""

import math
from typing import Sequence


def _mean(xs):
    return sum(xs) / len(xs)


def _std(xs, ddof=1):
    n = len(xs)
    if n - ddof <= 0:
        raise ValueError("need more observations than ddof")
    m = _mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (n - ddof))


def sharpe_ratio(returns: Sequence[float], risk_free=0.0,
                 periods_per_year=252) -> float:
    """Annualized Sharpe ratio of a periodic return series.

    ``(mean_excess / stdev) * sqrt(periods_per_year)`` where ``risk_free`` is the
    per-period risk-free return. Sample standard deviation (ddof=1).
    """
    if len(returns) < 2:
        raise ValueError("need at least two returns")
    excess = [r - risk_free for r in returns]
    sd = _std(excess)
    if sd <= 0.0:
        raise ValueError("zero-variance returns")
    return _mean(excess) / sd * math.sqrt(periods_per_year)


def sortino_ratio(returns: Sequence[float], risk_free=0.0, target=0.0,
                  periods_per_year=252) -> float:
    """Annualized Sortino ratio: excess mean over downside deviation.

    Downside deviation uses only returns below ``target`` (root-mean-square of
    the shortfalls, divided by the full sample count -- the standard
    convention). Raises if there is no downside.
    """
    if len(returns) < 2:
        raise ValueError("need at least two returns")
    excess_mean = _mean([r - risk_free for r in returns])
    downside = [min(r - target, 0.0) for r in returns]
    dd = math.sqrt(sum(d * d for d in downside) / len(returns))
    if dd <= 0.0:
        raise ValueError("no downside deviation (all returns >= target)")
    return excess_mean / dd * math.sqrt(periods_per_year)


def max_drawdown(returns: Sequence[float]) -> float:
    """Maximum peak-to-trough drawdown of the cumulative-return curve.

    Compounds the periodic returns into an equity curve and returns the largest
    fractional drop from a running peak, as a non-negative number (0.2 = a 20%
    drawdown). Empty or all-rising series give 0.
    """
    if not returns:
        return 0.0
    equity = 1.0
    peak = 1.0
    mdd = 0.0
    for r in returns:
        equity *= (1.0 + r)
        if equity > peak:
            peak = equity
        drop = (peak - equity) / peak
        if drop > mdd:
            mdd = drop
    return mdd


def drawdown_curve(returns: Sequence[float]) -> list:
    """Per-period underwater curve: fractional drop from the running peak.

    Compounds the returns into an equity curve and returns, for each period, the
    non-negative drawdown ``(peak - equity)/peak`` at that point (0 at a new
    high). The maximum of this curve is :func:`max_drawdown`.
    """
    out = []
    equity = 1.0
    peak = 1.0
    for r in returns:
        equity *= (1.0 + r)
        if equity > peak:
            peak = equity
        out.append((peak - equity) / peak)
    return out


def longest_drawdown_duration(returns: Sequence[float]) -> int:
    """Longest run of consecutive underwater periods (below a prior peak).

    Counts the maximum number of periods between a peak and the point the equity
    curve first recovers to (or exceeds) it. A series that never falls below its
    running peak returns 0.
    """
    equity = 1.0
    peak = 1.0
    longest = 0
    current = 0
    for r in returns:
        equity *= (1.0 + r)
        if equity >= peak:
            peak = equity
            current = 0
        else:
            current += 1
            if current > longest:
                longest = current
    return longest


def rolling_sharpe(returns: Sequence[float], window: int, risk_free=0.0,
                   periods_per_year=252) -> list:
    """Annualized Sharpe ratio over each trailing window of ``window`` periods.

    Returns one Sharpe per window position (``len(returns) - window + 1``
    values), each computed by :func:`sharpe_ratio` on that slice. A
    zero-variance window yields ``float('nan')`` rather than raising, so the
    series stays aligned.
    """
    n = len(returns)
    if window < 2:
        raise ValueError("window must be >= 2")
    if window > n:
        raise ValueError("window longer than the series")
    out = []
    for i in range(n - window + 1):
        chunk = returns[i:i + window]
        try:
            out.append(sharpe_ratio(chunk, risk_free, periods_per_year))
        except ValueError:
            out.append(float("nan"))
    return out


def calmar_ratio(returns: Sequence[float], periods_per_year=252) -> float:
    """Calmar ratio: annualized return divided by the maximum drawdown.

    Annualized return is the geometric ``(prod(1+r))^{periods_per_year/n} - 1``.
    Raises if there is no drawdown (undefined ratio).
    """
    if len(returns) < 2:
        raise ValueError("need at least two returns")
    growth = 1.0
    for r in returns:
        growth *= (1.0 + r)
    ann_return = growth ** (periods_per_year / len(returns)) - 1.0
    mdd = max_drawdown(returns)
    if mdd <= 0.0:
        raise ValueError("no drawdown (Calmar undefined)")
    return ann_return / mdd


def hit_rate(returns: Sequence[float]) -> float:
    """Fraction of periods with a strictly positive return."""
    if not returns:
        raise ValueError("need at least one return")
    return sum(1 for r in returns if r > 0.0) / len(returns)


def profit_factor(returns: Sequence[float]) -> float:
    """Gross profits divided by gross losses (absolute).

    Returns ``inf`` when there are no losing periods. Raises if there are no
    gains and no losses.
    """
    gains = sum(r for r in returns if r > 0.0)
    losses = -sum(r for r in returns if r < 0.0)
    if gains == 0.0 and losses == 0.0:
        raise ValueError("no gains or losses")
    if losses == 0.0:
        return float("inf")
    return gains / losses
