"""Portfolio rebalancing: weight drift, turnover, no-trade bands, cost drag.

Between rebalances, portfolio weights drift with the assets' returns. Rebalancing
back to a target incurs turnover and transaction costs. This module computes the
drifted weights, one-way turnover, the transaction-cost drag, and a no-trade-band
rebalance that only trades positions outside a tolerance. Pure standard library.
"""


def drift_weights(weights, asset_returns):
    """Buy-and-hold weights after one period of ``asset_returns``.

    Each position grows by ``(1 + r_i)``; the new weights are the grown values
    renormalized to sum to one. The starting point for the next rebalance decision.
    """
    n = len(weights)
    if len(asset_returns) != n:
        raise ValueError("weights and asset_returns must have equal length")
    grown = [weights[i] * (1.0 + asset_returns[i]) for i in range(n)]
    total = sum(grown)
    if total <= 0:
        raise ValueError("portfolio value must stay positive")
    return [g / total for g in grown]


def turnover(current_weights, target_weights):
    """One-way turnover ``0.5 * sum |target - current|`` (fraction of the book).

    The fraction of the portfolio traded to move from current to target weights;
    zero when already on target, up to one for a full turnover.
    """
    n = len(current_weights)
    if len(target_weights) != n:
        raise ValueError("weight vectors must have equal length")
    return 0.5 * sum(abs(target_weights[i] - current_weights[i]) for i in range(n))


def transaction_cost(current_weights, target_weights, cost_bps):
    """Transaction-cost drag of a rebalance: ``2 * turnover * cost_bps / 1e4``.

    Costs the round-trip (both sides) at ``cost_bps`` basis points of the traded
    notional. Zero when no trade is needed.
    """
    if cost_bps < 0:
        raise ValueError("cost_bps must be non-negative")
    return 2.0 * turnover(current_weights, target_weights) * cost_bps / 1e4


def no_trade_band_rebalance(current_weights, target_weights, band):
    """Rebalance only positions whose drift exceeds a ``band`` tolerance.

    Positions within ``band`` of their target are left untouched (no trade);
    those outside are moved to target. The remaining weight from the traded legs
    is left as-is (the untouched legs keep their drifted weight), so the result is
    renormalized to sum to one. Reduces turnover versus a full rebalance.
    """
    if band < 0:
        raise ValueError("band must be non-negative")
    n = len(current_weights)
    if len(target_weights) != n:
        raise ValueError("weight vectors must have equal length")
    new = []
    for i in range(n):
        if abs(current_weights[i] - target_weights[i]) > band:
            new.append(target_weights[i])
        else:
            new.append(current_weights[i])
    total = sum(new)
    if total <= 0:
        raise ValueError("weights must sum to a positive number")
    return [w / total for w in new]
