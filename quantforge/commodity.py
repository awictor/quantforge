"""Commodity forward pricing under the cost-of-carry model.

A storable commodity's forward is the spot compounded at the financing rate plus
storage cost and net of the convenience yield -- the benefit of holding the
physical good rather than a claim on it:

    F(T) = S * exp((r + u - y) * T)

with ``r`` the risk-free rate, ``u`` the (continuous, proportional) storage cost,
and ``y`` the convenience yield. When ``y > r + u`` the curve is in backwardation
(forwards below spot); otherwise it is in contango. This module prices the
forward, inverts the market forward for the implied convenience yield, builds a
forward curve, and classifies the curve shape. Pure standard library.
"""

import math


def commodity_forward(spot, r, maturity, storage_cost=0.0, convenience_yield=0.0):
    """Cost-of-carry forward ``S * exp((r + u - y) * T)``.

    ``storage_cost`` (u) and ``convenience_yield`` (y) are continuous proportional
    rates. Storage lifts the forward (a cost of carrying the physical), the
    convenience yield lowers it (a benefit of holding it).
    """
    if spot <= 0:
        raise ValueError("spot must be positive")
    if maturity < 0:
        raise ValueError("maturity must be non-negative")
    return spot * math.exp((r + storage_cost - convenience_yield) * maturity)


def implied_convenience_yield(spot, forward, r, maturity, storage_cost=0.0):
    """Convenience yield implied by a market forward (inverts the carry formula).

    Solves ``F = S exp((r + u - y) T)`` for ``y``:
    ``y = r + u - ln(F/S)/T``. Inverse of :func:`commodity_forward`.
    """
    if spot <= 0 or forward <= 0:
        raise ValueError("spot and forward must be positive")
    if maturity <= 0:
        raise ValueError("maturity must be positive")
    return r + storage_cost - math.log(forward / spot) / maturity


def implied_storage_cost(spot, forward, r, maturity, convenience_yield=0.0):
    """Storage cost implied by a market forward (inverts the carry formula).

    ``u = ln(F/S)/T - r + y``. Inverse of :func:`commodity_forward` in ``u``.
    """
    if spot <= 0 or forward <= 0:
        raise ValueError("spot and forward must be positive")
    if maturity <= 0:
        raise ValueError("maturity must be positive")
    return math.log(forward / spot) / maturity - r + convenience_yield


def net_cost_of_carry(r, storage_cost=0.0, convenience_yield=0.0):
    """Net proportional carry rate ``r + u - y`` (the forward's growth rate)."""
    return r + storage_cost - convenience_yield


def commodity_forward_curve(spot, r, maturities, storage_cost=0.0,
                            convenience_yield=0.0):
    """Forward prices across a list of maturities under one carry rate.

    Returns ``[(T, F(T)), ...]`` from :func:`commodity_forward` at each maturity.
    """
    return [(T, commodity_forward(spot, r, T, storage_cost, convenience_yield))
            for T in maturities]


def is_backwardation(r, storage_cost=0.0, convenience_yield=0.0):
    """True when the curve is in backwardation (``y > r + u``, forwards below spot).

    Backwardation occurs when the convenience yield exceeds the financing-plus-
    storage carry, so the net carry :func:`net_cost_of_carry` is negative and
    forwards fall with maturity.
    """
    return net_cost_of_carry(r, storage_cost, convenience_yield) < 0.0
