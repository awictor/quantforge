"""Realized kernel estimator of integrated variance (Barndorff-Nielsen et al.).

High-frequency realized variance is badly upward-biased by market-microstructure
noise (bid-ask bounce, discreteness). The realized kernel corrects this by adding
weighted autocovariances of the intraday returns,

    K = gamma_0 + sum_{h=1}^{H} k((h-1)/H) (gamma_h + gamma_{-h}),

where ``gamma_h`` is the ``h``-th return autocovariance and ``k(.)`` is a smooth
weight (the flat-top Parzen kernel here). The noise contribution in ``gamma_0`` is
cancelled by the weighted higher-lag terms, so the estimator is consistent for the
integrated variance even under noise. Pure standard library.
"""


def _parzen(x):
    """Parzen kernel weight ``k(x)`` on ``[0, 1]`` (``k(0)=1``, ``k(1)=0``)."""
    if x <= 0.5:
        return 1.0 - 6.0 * x * x + 6.0 * x ** 3
    return 2.0 * (1.0 - x) ** 3


def _autocov(returns, h):
    """``h``-th sample autocovariance of returns (not mean-adjusted; returns ~ 0)."""
    n = len(returns)
    return sum(returns[i] * returns[i - h] for i in range(h, n))


def realized_kernel(prices, bandwidth=None):
    """Realized-kernel integrated-variance estimate from a (log) price series.

    Uses the flat-top Parzen kernel over ``bandwidth`` lags of the intraday-return
    autocovariances. If ``bandwidth`` is ``None`` it defaults to the
    Barndorff-Nielsen rule of thumb ``H ~ n^{3/5}`` (capped below the number of
    returns). Robust to i.i.d. microstructure noise, unlike the naive realized
    variance, and always non-negative for the Parzen kernel. Requires at least three
    prices.
    """
    n_prices = len(prices)
    if n_prices < 3:
        raise ValueError("need at least 3 prices")
    returns = [prices[i + 1] - prices[i] for i in range(n_prices - 1)]
    n = len(returns)
    if bandwidth is None:
        bandwidth = max(1, int(round(n ** 0.6)))
    bandwidth = min(bandwidth, n - 1)

    total = _autocov(returns, 0)     # gamma_0
    for h in range(1, bandwidth + 1):
        weight = _parzen((h - 1) / bandwidth)
        gamma_h = _autocov(returns, h)
        total += weight * 2.0 * gamma_h     # gamma_h + gamma_{-h} = 2 gamma_h
    return max(0.0, total)
