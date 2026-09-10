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

    # Initial guess: Corrado-Miller (1996) rational approximation. It is
    # accurate away from the money too (not just ATM like Brenner-Subrahmanyam),
    # cutting Newton iterations several-fold. Falls back to Brenner-Subrahmanyam
    # if the discriminant goes negative (deep wings).
    disc = math.exp(-r * t)
    X = K * disc                       # present value of the strike
    # Express the quote as a call-equivalent for the Corrado-Miller formula.
    if ot is OptionType.CALL:
        c = target_price
    else:
        c = target_price + S * math.exp((b - r) * t) - X   # put-call parity
    S_disc = S * math.exp((b - r) * t)
    a = c - 0.5 * (S_disc - X)
    radicand = a * a - (S_disc - X) ** 2 / math.pi
    if radicand >= 0 and (S_disc + X) > 0:
        sigma = (math.sqrt(2.0 * math.pi / t) / (S_disc + X)) * (a + math.sqrt(radicand))
    else:
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


def implied_vol_smile(strikes, prices, S, t, r, option_type=OptionType.CALL,
                      b=None, forward=None):
    """Invert a whole chain of quotes to an implied-vol smile in one call.

    Args:
        strikes, prices: equal-length option-quote arrays at one expiry.
        forward: optional forward for the log-moneyness output; defaults to the
            carry-implied forward ``S e^{b t}``.

    Returns a list of ``(log_moneyness, implied_vol)`` pairs sorted by strike,
    skipping any quote outside the no-arbitrage band (those cannot be inverted).
    ``log_moneyness = ln(K / forward)``, the standard smile x-axis.
    """
    ot = _coerce_type(option_type)
    if b is None:
        b = r
    if len(strikes) != len(prices):
        raise ValueError("strikes and prices must be the same length")
    if forward is None:
        forward = S * math.exp(b * t)

    out = []
    for K, px in sorted(zip(strikes, prices)):
        try:
            iv = implied_volatility(px, S, K, t, r, ot, b)
        except ValueError:
            continue   # quote outside the no-arbitrage band; drop it
        out.append((math.log(K / forward), iv))
    return out
