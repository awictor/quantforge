"""Covered-call and cash-secured-put income analytics.

These are the standard yield/return metrics option-income traders quote, built
on the pricing engine so the premium is a model value (or pass a market premium
directly). All returns are simple (not compounded) and annualized by ``t``.
"""

import math
from dataclasses import dataclass

from .bsm import call_price, put_price


@dataclass(frozen=True)
class IncomeMetrics:
    premium: float           # option premium collected (per share)
    static_yield: float      # premium / capital, over the holding period
    annualized_yield: float  # static_yield / t
    if_assigned_return: float  # total simple return if assigned/exercised
    breakeven: float           # underlying level where the position breaks even


def covered_call(S, K, t, r, sigma, b=None, premium=None) -> IncomeMetrics:
    """Covered call: long stock at ``S``, short a call struck at ``K``.

    Args:
        premium: option premium; if None, uses the BSM call value.

    Yields are relative to the stock capital ``S``. If assigned (S_T >= K) the
    return is the capped gain to the strike plus the premium; the breakeven is
    ``S - premium`` (the stock can fall by the premium before a loss).
    """
    if b is None:
        b = r
    if premium is None:
        premium = call_price(S, K, t, r, sigma, b=b)
    static_yield = premium / S
    annualized = static_yield / t if t > 0 else 0.0
    # Assigned: sell stock at K, keep premium. Total P&L = (K - S) + premium.
    if_assigned = ((K - S) + premium) / S
    breakeven = S - premium
    return IncomeMetrics(premium=premium, static_yield=static_yield,
                         annualized_yield=annualized, if_assigned_return=if_assigned,
                         breakeven=breakeven)


def cash_secured_put(S, K, t, r, sigma, b=None, premium=None) -> IncomeMetrics:
    """Cash-secured put: short a put struck at ``K``, holding ``K`` cash.

    Args:
        premium: option premium; if None, uses the BSM put value.

    Yields are relative to the secured cash ``K``. If assigned (S_T <= K) the
    trader buys the stock at ``K`` net of premium, so the effective purchase and
    breakeven price is ``K - premium``.
    """
    if b is None:
        b = r
    if premium is None:
        premium = put_price(S, K, t, r, sigma, b=b)
    static_yield = premium / K
    annualized = static_yield / t if t > 0 else 0.0
    # If assigned the return over the secured cash is just the premium kept
    # (before the stock's own P&L); the effective buy price is K - premium.
    if_assigned = premium / K
    breakeven = K - premium
    return IncomeMetrics(premium=premium, static_yield=static_yield,
                         annualized_yield=annualized, if_assigned_return=if_assigned,
                         breakeven=breakeven)
