"""Order-flow microstructure measures estimated from trade data.

Complements the theoretical Kyle impact and the low-frequency liquidity proxies
with estimators built directly from signed trades:

- ``kyle_lambda_regression`` -- the price-impact coefficient from an OLS regression
  of price changes on signed order flow (the empirical counterpart of Kyle's lambda),
- ``order_flow_imbalance`` -- net signed volume over total volume in a window,
- ``vpin`` -- the volume-synchronized probability of informed trading (Easley,
  Lopez de Prado, O'Hara), the mean absolute order imbalance across equal-volume
  buckets.

Pure standard library.
"""


def kyle_lambda_regression(price_changes, signed_volumes):
    """Empirical Kyle's lambda: slope of price change on signed order flow.

    Fits ``dP = alpha + lambda * signed_volume`` by ordinary least squares and
    returns ``lambda`` (price move per unit of net signed volume). ``signed_volumes``
    should be signed by trade direction (buys positive, sells negative). A larger
    lambda means a less liquid, higher-impact market. Requires at least three
    aligned observations with non-constant order flow.
    """
    n = len(price_changes)
    if n < 3 or len(signed_volumes) != n:
        raise ValueError("need at least 3 aligned observations")
    mx = sum(signed_volumes) / n
    my = sum(price_changes) / n
    sxx = sum((signed_volumes[i] - mx) ** 2 for i in range(n))
    if sxx <= 0.0:
        raise ValueError("order flow has zero variance")
    sxy = sum((signed_volumes[i] - mx) * (price_changes[i] - my) for i in range(n))
    return sxy / sxx


def order_flow_imbalance(buy_volumes, sell_volumes):
    """Net signed volume over total volume: ``sum(buy - sell) / sum(buy + sell)``.

    In ``[-1, 1]``: positive when buys dominate, negative when sells do. Aligned
    non-negative volume series.
    """
    n = len(buy_volumes)
    if n == 0 or len(sell_volumes) != n:
        raise ValueError("buy and sell volumes must be equal-length, non-empty")
    if any(b < 0.0 for b in buy_volumes) or any(s < 0.0 for s in sell_volumes):
        raise ValueError("volumes must be non-negative")
    total = sum(buy_volumes) + sum(sell_volumes)
    if total <= 0.0:
        raise ValueError("total volume must be positive")
    return (sum(buy_volumes) - sum(sell_volumes)) / total


def vpin(buy_volumes, sell_volumes):
    """Volume-synchronized probability of informed trading (VPIN).

    Given per-bucket buy and sell volumes (equal-volume buckets), VPIN is the mean
    of ``|buy - sell| / (buy + sell)`` across buckets -- the average absolute order
    imbalance. In ``[0, 1]``: near zero when buys and sells balance, near one when
    trading is one-sided (a proxy for informed order flow / toxicity). Aligned
    non-negative series with at least one bucket.
    """
    n = len(buy_volumes)
    if n == 0 or len(sell_volumes) != n:
        raise ValueError("buy and sell volumes must be equal-length, non-empty")
    if any(b < 0.0 for b in buy_volumes) or any(s < 0.0 for s in sell_volumes):
        raise ValueError("volumes must be non-negative")
    total = 0.0
    count = 0
    for i in range(n):
        v = buy_volumes[i] + sell_volumes[i]
        if v > 0.0:
            total += abs(buy_volumes[i] - sell_volumes[i]) / v
            count += 1
    if count == 0:
        raise ValueError("all buckets have zero volume")
    return total / count
