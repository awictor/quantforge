"""Constant Proportion Portfolio Insurance (CPPI).

CPPI dynamically allocates between a risky asset and a safe asset to guarantee a
floor. The *cushion* is the wealth above the (discounted) floor; the risky
exposure is a fixed *multiplier* times the cushion, capped at total wealth and
floored at zero. As wealth approaches the floor the risky exposure shrinks to
zero, locking in the guarantee. This module computes the exposure, the discounted
floor, and simulates a wealth path. Pure standard library.
"""

import math


def discounted_floor(floor, r, horizon):
    """Present value of the guaranteed floor: ``floor * e^{-r * horizon}``.

    The bond floor CPPI must stay above today so the terminal wealth is at least
    ``floor``.
    """
    if horizon < 0:
        raise ValueError("horizon must be non-negative")
    return floor * math.exp(-r * horizon)


def cushion(wealth, floor_pv):
    """Cushion: wealth above the (discounted) floor, ``max(wealth - floor_pv, 0)``."""
    return max(wealth - floor_pv, 0.0)


def risky_exposure(wealth, floor_pv, multiplier):
    """CPPI risky-asset exposure ``clamp(multiplier * cushion, 0, wealth)``.

    The dollar amount in the risky asset: the multiplier times the
    :func:`cushion`, capped at total wealth (no leverage) and floored at zero.
    Zero once wealth hits the floor, protecting the guarantee.
    """
    if multiplier < 0:
        raise ValueError("multiplier must be non-negative")
    if wealth < 0:
        raise ValueError("wealth must be non-negative")
    e = multiplier * cushion(wealth, floor_pv)
    return min(max(e, 0.0), wealth)


def cppi_path(initial_wealth, floor, multiplier, risky_returns, r, dt):
    """Simulate a CPPI wealth path over a sequence of risky-asset returns.

    Each step: allocate :func:`risky_exposure` to the risky asset (rest at the safe
    rate ``r`` over ``dt``), apply that period's ``risky_return``, and roll forward.
    The floor is discounted to each step's remaining horizon. Returns the list of
    period-end wealths. Wealth stays at or above the floor for a multiplier within
    the gap-risk limit.
    """
    if initial_wealth <= 0:
        raise ValueError("initial_wealth must be positive")
    n = len(risky_returns)
    wealth = initial_wealth
    out = []
    safe_growth = math.exp(r * dt)
    for k, rr in enumerate(risky_returns):
        remaining = (n - k) * dt
        fpv = discounted_floor(floor, r, remaining)
        e = risky_exposure(wealth, fpv, multiplier)
        safe = wealth - e
        wealth = e * (1.0 + rr) + safe * safe_growth
        out.append(wealth)
    return out
