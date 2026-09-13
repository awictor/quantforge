"""Trade-sign classification: infer buyer/seller-initiated trades from tape data.

Order-flow measures (VPIN, order-flow imbalance, empirical Kyle's lambda) need each
trade signed as buyer-initiated (+1) or seller-initiated (-1), but raw tapes rarely
carry the aggressor side. The standard classifiers infer it:

- ``tick_rule`` -- sign by the last price change (uptick = buy, downtick = sell,
  carry the previous sign on a zero tick),
- ``quote_rule`` -- sign by the trade's side of the prevailing bid-ask midpoint,
- ``lee_ready`` -- the Lee-Ready (1991) hybrid: use the quote rule away from the
  midpoint and fall back to the tick rule at the midpoint.

Each returns a list of ``+1 / -1`` signs aligned with the trades. Pure standard
library.
"""


def tick_rule(prices):
    """Classify trades by the tick test on the trade-price series.

    ``+1`` on an uptick, ``-1`` on a downtick, and the previous sign carried forward
    on a zero tick (the first trade defaults to ``+1``). Needs at least one price.
    """
    if not prices:
        raise ValueError("need at least one price")
    signs = []
    last = 1
    prev_price = None
    for p in prices:
        if prev_price is None or p == prev_price:
            signs.append(last)
        else:
            last = 1 if p > prev_price else -1
            signs.append(last)
        prev_price = p
    return signs


def quote_rule(prices, bids, asks):
    """Classify trades by their side of the prevailing bid-ask midpoint.

    ``+1`` above the midpoint (buyer-initiated), ``-1`` below, and ``0`` exactly at
    the midpoint (unclassified; :func:`lee_ready` resolves these). Aligned series.
    """
    n = len(prices)
    if n == 0 or len(bids) != n or len(asks) != n:
        raise ValueError("prices, bids, asks must be equal-length, non-empty")
    signs = []
    for i in range(n):
        mid = 0.5 * (bids[i] + asks[i])
        if prices[i] > mid:
            signs.append(1)
        elif prices[i] < mid:
            signs.append(-1)
        else:
            signs.append(0)
    return signs


def lee_ready(prices, bids, asks):
    """Lee-Ready trade classification: quote rule with a tick-rule tiebreak.

    Trades away from the midpoint are signed by the quote rule; trades exactly at
    the midpoint are signed by the tick rule on the trade-price series. Returns a
    list of ``+1 / -1`` signs (no zeros). Aligned series of at least one trade.
    """
    n = len(prices)
    if n == 0 or len(bids) != n or len(asks) != n:
        raise ValueError("prices, bids, asks must be equal-length, non-empty")
    quotes = quote_rule(prices, bids, asks)
    ticks = tick_rule(prices)
    return [quotes[i] if quotes[i] != 0 else ticks[i] for i in range(n)]
