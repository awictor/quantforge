"""Digital-option super-replication by a vanilla call/put spread.

A cash-or-nothing digital paying ``Q`` if ``S_T > K`` (a call digital) has an
unbounded gamma/delta right at the strike near expiry, so it cannot be hedged
with the digital's own delta. Desks instead **over-hedge** it with a tight
vanilla call spread that *dominates* the digital payoff everywhere:

    long  Q/w calls struck at K - w
    short Q/w calls struck at K

Its payoff is 0 below ``K - w``, ramps to ``Q`` at ``K``, and stays ``Q`` above
- so it pays at least the digital's ``Q * 1{S>K}`` at every terminal spot. The
spread is what the seller actually trades; its cost is a conservative (upper-
bound) price for the digital, converging to the true digital value as the width
``w`` shrinks. The overhedge cost minus the fair digital value is the cushion
the desk keeps for the un-hedgeable pin risk.
"""

from dataclasses import dataclass

from .bsm import call_price, put_price
from .exotics import cash_or_nothing
from .bsm import OptionType, _coerce_type


@dataclass(frozen=True)
class Overhedge:
    cost: float             # price of the replicating spread (upper bound)
    digital_value: float    # fair value of the digital being hedged
    cushion: float          # cost - digital_value (the pin-risk buffer)
    long_strike: float
    short_strike: float
    quantity: float         # options per leg = Q / width


def digital_call_overhedge(S, K, t, r, sigma, cash=1.0, width=None, b=None):
    """Super-replicate a cash-or-nothing CALL digital with a call spread.

    Longs ``cash/width`` calls at ``K - width`` and shorts the same at ``K``, so
    the payoff dominates ``cash * 1{S_T > K}``. Returns an :class:`Overhedge`
    with the spread cost (a conservative price), the fair digital value, and the
    cushion between them.
    """
    if b is None:
        b = r
    if width is None:
        width = max(1e-3, 0.01 * K)   # default 1% strike width
    if width <= 0:
        raise ValueError("width must be positive")

    qty = cash / width
    long_k = K - width
    short_k = K
    long_leg = call_price(S, long_k, t, r, sigma, b=b)
    short_leg = call_price(S, short_k, t, r, sigma, b=b)
    cost = qty * (long_leg - short_leg)

    digital = cash_or_nothing(S, K, t, r, sigma, OptionType.CALL, b=b, cash=cash)
    return Overhedge(cost=cost, digital_value=digital, cushion=cost - digital,
                     long_strike=long_k, short_strike=short_k, quantity=qty)


def digital_put_overhedge(S, K, t, r, sigma, cash=1.0, width=None, b=None):
    """Super-replicate a cash-or-nothing PUT digital with a put spread.

    Longs ``cash/width`` puts at ``K + width`` and shorts the same at ``K``, so
    the payoff dominates ``cash * 1{S_T < K}``.
    """
    if b is None:
        b = r
    if width is None:
        width = max(1e-3, 0.01 * K)
    if width <= 0:
        raise ValueError("width must be positive")

    qty = cash / width
    long_k = K + width
    short_k = K
    long_leg = put_price(S, long_k, t, r, sigma, b=b)
    short_leg = put_price(S, short_k, t, r, sigma, b=b)
    cost = qty * (long_leg - short_leg)

    digital = cash_or_nothing(S, K, t, r, sigma, OptionType.PUT, b=b, cash=cash)
    return Overhedge(cost=cost, digital_value=digital, cushion=cost - digital,
                     long_strike=long_k, short_strike=short_k, quantity=qty)


def overhedge_payoff(oh: Overhedge, spot_at_expiry: float, is_call=True) -> float:
    """Terminal payoff of the replicating spread at ``spot_at_expiry``."""
    if is_call:
        lo, hi = oh.long_strike, oh.short_strike
        ramp = max(spot_at_expiry - lo, 0.0) - max(spot_at_expiry - hi, 0.0)
    else:
        # Put spread: long higher strike, short lower.
        hi, lo = oh.long_strike, oh.short_strike
        ramp = max(hi - spot_at_expiry, 0.0) - max(lo - spot_at_expiry, 0.0)
    return oh.quantity * ramp
