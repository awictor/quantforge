"""Leveraged / inverse ETF path and volatility drag.

A daily-rebalanced leveraged ETF delivers ``L`` times the *daily* return, not ``L``
times the period return. Compounding a constant leverage over volatile returns
creates a drag: the expected log return loses roughly ``0.5 L (L - 1) sigma^2`` per
period versus naive ``L`` times exposure. This module simulates the ETF path,
quantifies the drag, and gives the approximate multi-period return. Pure standard
library.
"""

import math


def leveraged_etf_path(underlying_returns, leverage, expense_ratio=0.0,
                       periods_per_year=252):
    """Daily-rebalanced leveraged ETF cumulative return path.

    Each period the ETF returns ``leverage * r - expense_ratio/periods_per_year``;
    the path compounds those. ``leverage`` may be negative (inverse ETFs). Returns
    the list of cumulative growth factors (starting after the first period).
    """
    daily_fee = expense_ratio / periods_per_year
    nav = 1.0
    out = []
    for r in underlying_returns:
        nav *= (1.0 + leverage * r - daily_fee)
        out.append(nav)
    return out


def volatility_drag(leverage, sigma):
    """Approximate per-period volatility drag of a leveraged ETF.

    ``0.5 * leverage * (leverage - 1) * sigma^2`` -- the expected log-return
    shortfall versus naive ``leverage`` times the underlying's log return, from the
    daily-rebalancing compounding. Zero at ``leverage`` 0 or 1; positive (a drag)
    for ``leverage > 1`` or ``leverage < 0``.
    """
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    return 0.5 * leverage * (leverage - 1.0) * sigma * sigma


def expected_leveraged_return(underlying_return, leverage, sigma, periods):
    """Approximate multi-period leveraged-ETF log return with drag.

    ``periods * (leverage * mu - drag)`` where ``mu`` is the per-period underlying
    log return and ``drag`` is :func:`volatility_drag`. Below naive
    ``periods * leverage * mu`` whenever there is drag.
    """
    if periods < 0:
        raise ValueError("periods must be non-negative")
    drag = volatility_drag(leverage, sigma)
    return periods * (leverage * underlying_return - drag)


def flat_market_decay(leverage, returns):
    """Cumulative leveraged return over a *round-trip* (net-flat) return path.

    A market that ends where it started but moved in between: the leveraged ETF
    still loses value to the drag. Returns the final growth factor minus one; below
    zero for ``|leverage| > 1`` over a volatile flat path.
    """
    path = leveraged_etf_path(returns, leverage)
    return path[-1] - 1.0 if path else 0.0
