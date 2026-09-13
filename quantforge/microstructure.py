"""Order-flow microstructure measures estimated from trade data.

Complements the theoretical Kyle impact and the low-frequency liquidity proxies
with estimators built directly from signed trades:

- ``kyle_lambda_regression`` -- the price-impact coefficient from an OLS regression
  of price changes on signed order flow (the empirical counterpart of Kyle's lambda),
- ``order_flow_imbalance`` -- net signed volume over total volume in a window,
- ``vpin`` -- the volume-synchronized probability of informed trading (Easley,
  Lopez de Prado, O'Hara), the mean absolute order imbalance across equal-volume
  buckets,
- the transaction-cost spread decomposition ``quoted_spread``,
  ``effective_spread``, ``realized_spread`` and ``price_impact``, which satisfy the
  identity ``effective = realized + price_impact``.

Pure standard library.
"""

import math


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


def _spread_inputs(trade_prices, mids, signs):
    n = len(trade_prices)
    if n == 0 or len(mids) != n or len(signs) != n:
        raise ValueError("trade_prices, mids, signs must be equal-length, non-empty")
    if any(s not in (-1, 1) for s in signs):
        raise ValueError("signs must be +1 (buy) or -1 (sell)")
    return n


def effective_spread(trade_prices, mids, signs):
    """Average effective (proportional) spread ``2 * sign * (price - mid) / mid``.

    The cost actually paid relative to the midpoint at the time of the trade:
    ``sign`` is ``+1`` for buys and ``-1`` for sells. Returned as the mean over the
    trades, in the same units as ``price / mid`` (a fraction). Wider than the quoted
    spread when trades walk the book, tighter when they occur inside it.
    """
    n = _spread_inputs(trade_prices, mids, signs)
    if any(m <= 0.0 for m in mids):
        raise ValueError("midpoints must be positive")
    return sum(2.0 * signs[i] * (trade_prices[i] - mids[i]) / mids[i]
               for i in range(n)) / n


def realized_spread(trade_prices, mids, future_mids, signs):
    """Average realized (proportional) spread ``2 * sign * (price - mid_future) / mid``.

    The portion of the effective spread the liquidity provider *keeps* -- the trade
    price against the midpoint a short horizon later (``future_mids``), so it nets
    out the permanent price move. ``mids`` is the quote midpoint at the trade,
    ``future_mids`` the midpoint after the impact horizon. Mean over the trades.
    """
    n = _spread_inputs(trade_prices, mids, signs)
    if len(future_mids) != n:
        raise ValueError("future_mids must align with the trades")
    if any(m <= 0.0 for m in mids):
        raise ValueError("midpoints must be positive")
    return sum(2.0 * signs[i] * (trade_prices[i] - future_mids[i]) / mids[i]
               for i in range(n)) / n


def price_impact(trade_prices, mids, future_mids, signs):
    """Average (proportional) price impact ``2 * sign * (mid_future - mid) / mid``.

    The permanent midpoint move in the trade's direction over the impact horizon --
    the informational half of the spread. By construction
    ``effective = realized + price_impact`` term by term (both defined with the same
    ``2 * sign / mid`` scaling), so this equals the gap between the effective and
    realized spreads.
    """
    n = _spread_inputs(trade_prices, mids, signs)
    if len(future_mids) != n:
        raise ValueError("future_mids must align with the trades")
    if any(m <= 0.0 for m in mids):
        raise ValueError("midpoints must be positive")
    return sum(2.0 * signs[i] * (future_mids[i] - mids[i]) / mids[i]
               for i in range(n)) / n


def quoted_spread(bids, asks):
    """Average proportional quoted spread ``(ask - bid) / midpoint``.

    The posted cost of a round trip, independent of where trades actually print.
    Aligned bid/ask series with positive midpoints.
    """
    n = len(bids)
    if n == 0 or len(asks) != n:
        raise ValueError("bids and asks must be equal-length, non-empty")
    total = 0.0
    for i in range(n):
        mid = 0.5 * (bids[i] + asks[i])
        if mid <= 0.0:
            raise ValueError("midpoints must be positive")
        total += (asks[i] - bids[i]) / mid
    return total / n
