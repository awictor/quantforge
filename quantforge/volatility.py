"""Realized (historical) volatility estimators from a price series.

All estimators return an *annualized* volatility. The range-based estimators
(Parkinson, Garman-Klass, Rogers-Satchell, Yang-Zhang) use intraday
high/low/open/close bars and are far more efficient than close-to-close for the
same sample size; Yang-Zhang additionally handles overnight jumps and drift.

Inputs are plain sequences of floats — no NumPy. ``periods_per_year`` defaults
to 252 (trading days); pass 252*390 for minute bars, etc.

References: Parkinson (1980), Garman-Klass (1980), Rogers-Satchell (1991),
Yang-Zhang (2000).
"""

import math
from dataclasses import dataclass
from typing import Sequence

_LN2 = math.log(2.0)


def _log_returns(closes: Sequence[float]):
    out = []
    for i in range(1, len(closes)):
        if closes[i] <= 0 or closes[i - 1] <= 0:
            raise ValueError("prices must be positive")
        out.append(math.log(closes[i] / closes[i - 1]))
    return out


def close_to_close(closes: Sequence[float], periods_per_year: int = 252,
                   ddof: int = 1) -> float:
    """Classic close-to-close realized volatility (annualized).

    Uses the sample standard deviation of log returns with ``ddof`` degrees of
    freedom removed (1 = unbiased sample variance).
    """
    rets = _log_returns(closes)
    n = len(rets)
    if n <= ddof:
        raise ValueError("not enough returns for the requested ddof")
    mean = sum(rets) / n
    var = sum((r - mean) ** 2 for r in rets) / (n - ddof)
    return math.sqrt(var * periods_per_year)


def ewma_vol(closes: Sequence[float], lam: float = 0.94,
             periods_per_year: int = 252) -> float:
    """RiskMetrics-style exponentially weighted volatility (annualized).

    Variance_t = lam * Variance_{t-1} + (1 - lam) * r_t^2, seeded with the
    first squared return. ``lam=0.94`` is the RiskMetrics daily default.
    """
    if not (0.0 < lam < 1.0):
        raise ValueError("lam must be in (0, 1)")
    rets = _log_returns(closes)
    if not rets:
        raise ValueError("need at least two prices")
    var = rets[0] ** 2
    for r in rets[1:]:
        var = lam * var + (1.0 - lam) * r * r
    return math.sqrt(var * periods_per_year)


def _validate_ohlc(opens, highs, lows, closes):
    n = len(closes)
    if not (len(opens) == len(highs) == len(lows) == n):
        raise ValueError("OHLC series must be the same length")
    if n < 1:
        raise ValueError("need at least one bar")
    for o, h, l, c in zip(opens, highs, lows, closes):
        if min(o, h, l, c) <= 0:
            raise ValueError("prices must be positive")
        if h < max(o, c) or l > min(o, c):
            raise ValueError("bar violates high/low bounds")
    return n


def parkinson(highs: Sequence[float], lows: Sequence[float],
              periods_per_year: int = 252) -> float:
    """Parkinson high-low range estimator (annualized).

    var = (1 / (4 ln2)) * mean( ln(H/L)^2 ). ~5x more efficient than
    close-to-close but ignores drift and overnight moves.
    """
    n = len(highs)
    if n < 1 or len(lows) != n:
        raise ValueError("need matching high/low series")
    acc = 0.0
    for h, l in zip(highs, lows):
        if h <= 0 or l <= 0 or h < l:
            raise ValueError("invalid high/low bar")
        acc += math.log(h / l) ** 2
    var = acc / (4.0 * _LN2 * n)
    return math.sqrt(var * periods_per_year)


def garman_klass(opens, highs, lows, closes, periods_per_year: int = 252) -> float:
    """Garman-Klass OHLC estimator (annualized).

    var = mean( 0.5*ln(H/L)^2 - (2ln2 - 1)*ln(C/O)^2 ). Uses the full bar; more
    efficient than Parkinson, still assumes no overnight jump or drift.
    """
    n = _validate_ohlc(opens, highs, lows, closes)
    acc = 0.0
    for o, h, l, c in zip(opens, highs, lows, closes):
        hl = math.log(h / l)
        co = math.log(c / o)
        acc += 0.5 * hl * hl - (2.0 * _LN2 - 1.0) * co * co
    var = acc / n
    return math.sqrt(var * periods_per_year)


def rogers_satchell(opens, highs, lows, closes, periods_per_year: int = 252) -> float:
    """Rogers-Satchell OHLC estimator (annualized).

    var = mean( ln(H/C)ln(H/O) + ln(L/C)ln(L/O) ). Drift-independent: stays
    unbiased even when the underlying has a non-zero mean return.
    """
    n = _validate_ohlc(opens, highs, lows, closes)
    acc = 0.0
    for o, h, l, c in zip(opens, highs, lows, closes):
        acc += (math.log(h / c) * math.log(h / o)
                + math.log(l / c) * math.log(l / o))
    var = acc / n
    return math.sqrt(var * periods_per_year)


def yang_zhang(opens, highs, lows, closes, periods_per_year: int = 252) -> float:
    """Yang-Zhang estimator (annualized): drift-independent and jump-robust.

    Combines overnight (close-to-open) variance, open-to-close variance, and
    the Rogers-Satchell term:

        var = var_overnight + k * var_open_to_close + (1 - k) * var_RS
        k   = 0.34 / (1.34 + (N+1)/(N-1))

    Requires the previous close, so bars are chained: overnight return uses
    ln(O_t / C_{t-1}).
    """
    n = _validate_ohlc(opens, highs, lows, closes)
    if n < 2:
        raise ValueError("Yang-Zhang needs at least two bars")

    # Overnight (close-to-open) log returns, chained across bars.
    overnight = [math.log(opens[i] / closes[i - 1]) for i in range(1, n)]
    # Open-to-close log returns on the same bars.
    open_close = [math.log(closes[i] / opens[i]) for i in range(1, n)]
    m = len(overnight)  # = n - 1

    mean_on = sum(overnight) / m
    var_on = sum((x - mean_on) ** 2 for x in overnight) / (m - 1) if m > 1 else 0.0
    mean_oc = sum(open_close) / m
    var_oc = sum((x - mean_oc) ** 2 for x in open_close) / (m - 1) if m > 1 else 0.0

    # Rogers-Satchell over the same bars (indices 1..n-1).
    rs = 0.0
    for i in range(1, n):
        o, h, l, c = opens[i], highs[i], lows[i], closes[i]
        rs += (math.log(h / c) * math.log(h / o)
               + math.log(l / c) * math.log(l / o))
    var_rs = rs / m

    k = 0.34 / (1.34 + (m + 1) / (m - 1)) if m > 1 else 0.34
    var = var_on + k * var_oc + (1.0 - k) * var_rs
    return math.sqrt(var * periods_per_year)


@dataclass(frozen=True)
class VolReport:
    close_to_close: float
    parkinson: float
    garman_klass: float
    rogers_satchell: float
    yang_zhang: float
    ewma: float


def vol_report(opens, highs, lows, closes, periods_per_year: int = 252,
               ewma_lambda: float = 0.94) -> VolReport:
    """Compute every estimator at once for an OHLC series."""
    return VolReport(
        close_to_close=close_to_close(closes, periods_per_year),
        parkinson=parkinson(highs, lows, periods_per_year),
        garman_klass=garman_klass(opens, highs, lows, closes, periods_per_year),
        rogers_satchell=rogers_satchell(opens, highs, lows, closes, periods_per_year),
        yang_zhang=yang_zhang(opens, highs, lows, closes, periods_per_year),
        ewma=ewma_vol(closes, ewma_lambda, periods_per_year),
    )


@dataclass(frozen=True)
class VolConePoint:
    window: int              # rolling window length (in bars)
    minimum: float
    p25: float
    median: float
    p75: float
    maximum: float
    current: float           # most-recent window's realized vol


def _percentile(sorted_vals, p):
    """Linear-interpolated percentile of an already-sorted list (p in [0,1])."""
    n = len(sorted_vals)
    if n == 1:
        return sorted_vals[0]
    idx = p * (n - 1)
    lo = int(math.floor(idx))
    hi = min(lo + 1, n - 1)
    frac = idx - lo
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[hi] * frac


def vol_cone(closes: Sequence[float], windows: Sequence[int],
             periods_per_year: int = 252):
    """Realized-volatility cone: the distribution of rolling realized vol per window.

    For each window length, computes the annualized close-to-close realized vol
    over every rolling block of returns of that length, then reports the min,
    25th/50th/75th percentiles, max, and the most-recent (current) value. This
    is the standard "vol cone" used to judge whether current realized vol is
    high or low versus its own history at each horizon.

    Args:
        closes: the price series.
        windows: rolling window lengths in *returns* (e.g. [5, 21, 63, 126]).
        periods_per_year: annualization factor.

    Returns a list of :class:`VolConePoint`, one per window (skipping windows
    too long for the data).
    """
    rets = _log_returns(closes)
    n = len(rets)
    out = []
    for w in sorted(windows):
        if w < 2 or w > n:
            continue
        vols = []
        for start in range(0, n - w + 1):
            block = rets[start:start + w]
            mean = sum(block) / w
            var = sum((x - mean) ** 2 for x in block) / (w - 1)
            vols.append(math.sqrt(var * periods_per_year))
        vols_sorted = sorted(vols)
        out.append(VolConePoint(
            window=w,
            minimum=vols_sorted[0],
            p25=_percentile(vols_sorted, 0.25),
            median=_percentile(vols_sorted, 0.5),
            p75=_percentile(vols_sorted, 0.75),
            maximum=vols_sorted[-1],
            current=vols[-1],   # most-recent rolling block
        ))
    return out
