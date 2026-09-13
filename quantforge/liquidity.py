"""Empirical liquidity and transaction-cost measures from price/volume data.

Low-frequency proxies that recover trading costs and price impact without an order
book:

- ``roll_spread`` -- Roll's (1984) effective spread from the serial covariance of
  price changes (bid-ask bounce induces negative first-order autocovariance),
- ``amihud_illiquidity`` -- Amihud's (2002) ratio of absolute return to dollar
  volume, the daily price impact per dollar traded,
- ``corwin_schultz_spread`` -- the Corwin-Schultz (2012) high-low spread estimator,
  which backs the spread out of two-day high/low ranges.

Pure standard library.
"""

import math


def roll_spread(prices):
    """Roll's implied effective spread from a series of transaction prices.

    Bid-ask bounce makes successive price changes negatively autocovaried; Roll's
    estimator is ``spread = 2 sqrt(-cov(dP_t, dP_{t-1}))``. When the sample
    autocovariance is non-negative (no detectable bounce) the estimate is zero.
    ``prices`` are levels; at least three are needed.
    """
    n = len(prices)
    if n < 3:
        raise ValueError("need at least 3 prices")
    dp = [prices[i + 1] - prices[i] for i in range(n - 1)]
    m = len(dp)
    mean = sum(dp) / m
    cov = sum((dp[i] - mean) * (dp[i - 1] - mean) for i in range(1, m)) / (m - 1)
    if cov >= 0.0:
        return 0.0
    return 2.0 * math.sqrt(-cov)


def amihud_illiquidity(returns, dollar_volumes):
    """Amihud illiquidity: average of ``|return| / dollar_volume`` over the period.

    A larger value means a given dollar of trading moves the price more (a less
    liquid asset). ``returns`` and ``dollar_volumes`` are aligned daily series.
    """
    n = len(returns)
    if n == 0 or len(dollar_volumes) != n:
        raise ValueError("returns and dollar_volumes must be equal-length, non-empty")
    if any(v <= 0.0 for v in dollar_volumes):
        raise ValueError("dollar volumes must be positive")
    return sum(abs(returns[i]) / dollar_volumes[i] for i in range(n)) / n


def corwin_schultz_spread(highs, lows):
    """Corwin-Schultz high-low bid-ask spread estimator.

    Uses consecutive daily high/low ranges: the two-day range reflects both
    volatility and the spread, while single-day ranges scale with volatility alone,
    so their combination isolates the spread. Returns the average estimated
    proportional spread over the sample, floored at zero each day (negative daily
    estimates, a known small-sample artefact, are set to zero). ``highs`` and
    ``lows`` are aligned and at least two long.
    """
    n = len(highs)
    if n < 2 or len(lows) != n:
        raise ValueError("highs and lows must be equal-length with at least 2 days")
    if any(highs[i] < lows[i] for i in range(n)):
        raise ValueError("each high must be at least its low")
    if any(lows[i] <= 0.0 for i in range(n)):
        raise ValueError("prices must be positive")

    k = 3.0 - 2.0 * math.sqrt(2.0)
    spreads = []
    for t in range(1, n):
        beta = (math.log(highs[t - 1] / lows[t - 1]) ** 2
                + math.log(highs[t] / lows[t]) ** 2)
        hi2 = max(highs[t - 1], highs[t])
        lo2 = min(lows[t - 1], lows[t])
        gamma = math.log(hi2 / lo2) ** 2
        alpha = ((math.sqrt(2.0 * beta) - math.sqrt(beta)) / k
                 - math.sqrt(gamma / k))
        s = 2.0 * (math.exp(alpha) - 1.0) / (1.0 + math.exp(alpha))
        spreads.append(max(0.0, s))
    return sum(spreads) / len(spreads)
