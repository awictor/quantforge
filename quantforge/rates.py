"""Interest-rate caps, floors, and collars priced on the Bachelier (normal) model.

A **caplet** pays ``tau * max(F - K, 0)`` at the end of an accrual period, where
``F`` is the forward rate over that period, ``K`` the cap strike, and ``tau`` the
year fraction (accrual). It is a call on the forward rate, discounted by the
period's zero-coupon bond. A **cap** is a strip of caplets across successive
periods; a **floor** is the strip of floorlets (puts); a **collar** is long a
cap and short a floor.

Rates can be negative, so the caplets use the Bachelier (normal-vol) model
rather than lognormal Black-76. Each period is described by its forward rate,
its normal volatility, its accrual, its expiry (when the rate sets), and its
discount factor to the payment date.
"""

import math
from dataclasses import dataclass
from typing import Sequence, List

from .bachelier import bachelier_price
from .bsm import OptionType


@dataclass(frozen=True)
class CapletPeriod:
    forward: float          # forward rate for the period
    expiry: float           # time (years) until the rate is observed/set
    accrual: float          # year fraction tau of the accrual period
    discount: float         # discount factor to the payment date
    sigma_n: float          # normal (absolute) volatility of the forward rate


def caplet_price(period: CapletPeriod, strike: float, is_cap: bool = True) -> float:
    """Price a single caplet (cap) or floorlet (floor).

    Value = discount * accrual * Bachelier(F, K, expiry, r=0, sigma_n),
    with the option being a call for a caplet and a put for a floorlet. The
    Bachelier price is taken undiscounted (r=0) and discounted explicitly by the
    period's bond factor, which is the market convention.
    """
    ot = OptionType.CALL if is_cap else OptionType.PUT
    if period.expiry <= 0:
        intrinsic = (max(period.forward - strike, 0.0) if is_cap
                     else max(strike - period.forward, 0.0))
        return period.discount * period.accrual * intrinsic
    undiscounted = bachelier_price(period.forward, strike, period.expiry, 0.0,
                                   period.sigma_n, ot)
    return period.discount * period.accrual * undiscounted


def cap_price(periods: Sequence[CapletPeriod], strike: float) -> float:
    """Price an interest-rate cap as the sum of its caplets."""
    return sum(caplet_price(p, strike, is_cap=True) for p in periods)


def floor_price(periods: Sequence[CapletPeriod], strike: float) -> float:
    """Price an interest-rate floor as the sum of its floorlets."""
    return sum(caplet_price(p, strike, is_cap=False) for p in periods)


def collar_price(periods: Sequence[CapletPeriod], cap_strike: float,
                 floor_strike: float) -> float:
    """Price a collar: long a cap at ``cap_strike``, short a floor at ``floor_strike``.

    The net value is ``cap - floor``; a zero-cost collar is the pair of strikes
    that makes this zero.
    """
    return cap_price(periods, cap_strike) - floor_price(periods, floor_strike)


def caplet_floorlet_parity(period: CapletPeriod, strike: float) -> float:
    """Caplet - floorlet at the same strike = discounted forward-minus-strike.

    A put-call-parity identity used to check the pricer:
    ``caplet - floorlet = discount * accrual * (F - K)``.
    """
    return period.discount * period.accrual * (period.forward - strike)
