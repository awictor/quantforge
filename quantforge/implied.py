"""Implied-volatility solver.

Newton-Raphson on vega with a Brent/bisection fallback so the solver is both
fast for well-behaved inputs and robust in the wings where vega collapses.
"""

import math

from .bsm import price, vega, OptionType, _coerce_type


def _no_arbitrage_bounds(S, K, t, r, b, ot):
    """European price bounds; used to reject arbitrage-violating quotes."""
    carry = math.exp((b - r) * t)
    disc = math.exp(-r * t)
    fwd_disc = S * carry
    if ot is OptionType.CALL:
        lower = max(fwd_disc - K * disc, 0.0)
        upper = fwd_disc
    else:
        lower = max(K * disc - fwd_disc, 0.0)
        upper = K * disc
    return lower, upper


def implied_volatility(
    target_price, S, K, t, r, option_type=OptionType.CALL, b=None,
    tol=1e-8, max_iter=100, lo=1e-9, hi=10.0,
):
    """Solve for the volatility that reproduces ``target_price``.

    Returns the implied vol, or raises ValueError if the quote is outside the
    no-arbitrage band (no finite vol can produce it).
    """
    ot = _coerce_type(option_type)
    if b is None:
        b = r
    if t <= 0:
        raise ValueError("cannot imply vol at or past expiry")

    lower, upper = _no_arbitrage_bounds(S, K, t, r, b, ot)
    # Allow a tiny epsilon for floating-point slop at the boundaries.
    eps = 1e-12
    if target_price < lower - eps or target_price > upper + eps:
        raise ValueError(
            f"price {target_price} outside no-arbitrage band "
            f"[{lower:.6g}, {upper:.6g}]"
        )
    if target_price <= lower + eps:
        return lo
    if target_price >= upper - eps:
        return hi

    # Initial guess: Brenner-Subrahmanyam ATM approximation.
    sigma = math.sqrt(2.0 * math.pi / t) * target_price / S
    sigma = min(max(sigma, lo), hi)

    f_lo = price(S, K, t, r, lo, ot, b) - target_price
    f_hi = price(S, K, t, r, hi, ot, b) - target_price
    if f_lo * f_hi > 0:
        # Should not happen given the band check, but guard anyway.
        raise ValueError("failed to bracket implied volatility")

    a, fa = lo, f_lo
    bb, fb = hi, f_hi
    for _ in range(max_iter):
        v = price(S, K, t, r, sigma, ot, b) - target_price
        if abs(v) < tol:
            return sigma
        # Maintain the bracket.
        if fa * v < 0:
            bb, fb = sigma, v
        else:
            a, fa = sigma, v

        vg = vega(S, K, t, r, sigma, b)
        step_ok = False
        if vg > 1e-12:
            newton = sigma - v / vg
            if a < newton < bb:
                sigma = newton
                step_ok = True
        if not step_ok:
            sigma = 0.5 * (a + bb)  # bisection fallback

    return sigma
