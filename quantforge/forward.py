"""Implied forward and dividend extraction from an option chain.

Put-call parity for European options states, at each strike ``K``,

    C(K) - P(K) = D * (F - K)

where ``D = exp(-r t)`` is the discount factor and ``F`` the forward price of
the underlying. Given a chain of call and put mid prices across strikes at one
expiry, the points ``(K, C - K)`` versus ``P`` are linear, so a least-squares
line recovers both the forward ``F`` and the discount factor ``D`` at once --
no volatility assumption needed. From the forward we back out the market's
implied dividend yield (or borrow cost).

This is the standard desk technique for building an arbitrage-consistent
forward curve directly from listed option prices.
"""

import math
from dataclasses import dataclass
from typing import Sequence


@dataclass(frozen=True)
class ForwardResult:
    forward: float          # implied forward price F
    discount_factor: float  # implied D = exp(-r t)
    implied_rate: float     # continuously-compounded r implied by D
    implied_div_yield: float  # continuous q such that F = S * exp((r - q) t)
    n_strikes: int
    rmse: float             # fit residual in price units


def implied_forward(strikes: Sequence[float], calls: Sequence[float],
                    puts: Sequence[float], t: float, spot: float = None):
    """Extract the implied forward and discount factor from a parity fit.

    Solves ``C - P = D*F - D*K`` as a straight line in ``K`` by ordinary least
    squares: the slope is ``-D`` and the intercept is ``D*F``.

    Args:
        strikes, calls, puts: equal-length chains at a single expiry.
        t: time to expiry in years (used to annualize the implied rate).
        spot: if given, also returns the implied continuous dividend yield.

    Returns a :class:`ForwardResult`.
    """
    n = len(strikes)
    if not (len(calls) == len(puts) == n):
        raise ValueError("strikes, calls and puts must be the same length")
    if n < 2:
        raise ValueError("need at least two strikes to fit a line")
    if t <= 0:
        raise ValueError("t must be positive")

    y = [calls[i] - puts[i] for i in range(n)]   # C - P
    x = list(strikes)

    # Ordinary least squares for y = a + m*x, where m = -D, a = D*F.
    xbar = sum(x) / n
    ybar = sum(y) / n
    sxx = sum((xi - xbar) ** 2 for xi in x)
    sxy = sum((x[i] - xbar) * (y[i] - ybar) for i in range(n))
    if sxx == 0:
        raise ValueError("all strikes identical; cannot fit")
    m = sxy / sxx
    a = ybar - m * xbar

    D = -m
    if D <= 0:
        raise ValueError("fit gave a non-positive discount factor; check inputs")
    F = a / D

    # Residuals in price units.
    ss = 0.0
    for i in range(n):
        pred = a + m * x[i]
        ss += (y[i] - pred) ** 2
    rmse = math.sqrt(ss / n)

    implied_rate = -math.log(D) / t

    div_yield = float("nan")
    if spot is not None:
        if spot <= 0:
            raise ValueError("spot must be positive")
        # F = S * exp((r - q) t)  ->  q = r - ln(F/S)/t
        div_yield = implied_rate - math.log(F / spot) / t

    return ForwardResult(forward=F, discount_factor=D, implied_rate=implied_rate,
                         implied_div_yield=div_yield, n_strikes=n, rmse=rmse)


def dividend_curve(chain_by_expiry, spot):
    """Bootstrap an implied dividend-yield term structure from a multi-expiry chain.

    Args:
        chain_by_expiry: iterable of ``(t, strikes, calls, puts)`` tuples, one
            per expiry.
        spot: current underlying spot.

    Returns a list of ``(t, ForwardResult)`` pairs sorted by expiry, each from
    :func:`implied_forward`. The ``implied_div_yield`` field of each result is
    the continuous dividend yield to that expiry (a point on the dividend curve).
    """
    out = []
    for t, strikes, calls, puts in chain_by_expiry:
        res = implied_forward(strikes, calls, puts, t, spot=spot)
        out.append((t, res))
    out.sort(key=lambda x: x[0])
    return out
