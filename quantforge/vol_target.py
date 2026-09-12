"""Volatility targeting: scale exposure to a constant risk level.

A volatility-targeting overlay scales position size inversely with recent
realized volatility so the portfolio runs at a roughly constant target vol:
leverage ``= target_vol / realized_vol``, capped at a maximum. This dampens risk
in turbulent regimes and levers up in calm ones. This module computes the
leverage and applies a rolling overlay to a return series. Pure standard library.
"""

import math


def target_leverage(target_vol, realized_vol, max_leverage=None):
    """Volatility-target leverage ``target_vol / realized_vol`` (capped).

    Above one when realized vol is below target (lever up), below one when it is
    above (de-risk). Capped at ``max_leverage`` if given, and floored at zero.
    """
    if target_vol < 0:
        raise ValueError("target_vol must be non-negative")
    if realized_vol <= 0:
        raise ValueError("realized_vol must be positive")
    lev = target_vol / realized_vol
    if max_leverage is not None:
        lev = min(lev, max_leverage)
    return max(lev, 0.0)


def _realized_vol(window, periods_per_year):
    n = len(window)
    if n < 2:
        return 0.0
    m = sum(window) / n
    var = sum((x - m) ** 2 for x in window) / (n - 1)
    return math.sqrt(var * periods_per_year)


def vol_targeted_returns(returns, target_vol, lookback, periods_per_year=252,
                         max_leverage=None):
    """Apply a rolling volatility-targeting overlay to a return series.

    For each period past the first ``lookback``, sizes the position at
    :func:`target_leverage` using the trailing ``lookback``-window annualized
    realized vol, and scales that period's return. Returns the overlaid return
    series (length ``len(returns) - lookback``). The overlay's realized vol sits
    near ``target_vol`` when the estimate tracks the true vol.
    """
    n = len(returns)
    if lookback < 2 or lookback >= n:
        raise ValueError("lookback must be in [2, len(returns) - 1]")
    if target_vol < 0:
        raise ValueError("target_vol must be non-negative")
    out = []
    for t in range(lookback, n):
        window = returns[t - lookback:t]
        rv = _realized_vol(window, periods_per_year)
        if rv <= 0.0:
            out.append(0.0)
            continue
        lev = target_leverage(target_vol, rv, max_leverage)
        out.append(lev * returns[t])
    return out


def realized_annualized_vol(returns, periods_per_year=252):
    """Annualized realized volatility of a return series (sample std)."""
    return _realized_vol(list(returns), periods_per_year)
