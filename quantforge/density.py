"""Breeden-Litzenberger risk-neutral density extraction.

Breeden and Litzenberger (1978) showed the risk-neutral probability density of
the underlying at expiry is recovered from the second derivative of the
(undiscounted) call price with respect to strike:

    f(K) = e^{r t} * d^2 C / dK^2

and the risk-neutral CDF from the first derivative:

    F(K) = 1 + e^{r t} * dC/dK.

Given a discrete call-price curve across strikes, we approximate these
derivatives by finite differences. The resulting density integrates to ~1 and,
combined with the discounted expectation identity, reprices any European payoff
by quadrature. This module exposes the density/CDF on a strike grid and a
generic payoff pricer built on it.
"""

import math
from typing import Sequence, Callable, List, Tuple


def risk_neutral_density(strikes: Sequence[float], calls: Sequence[float],
                         t: float, r: float) -> Tuple[List[float], List[float]]:
    """Estimate the risk-neutral pdf at the interior strikes.

    Uses a non-uniform central second difference of the call curve, so strikes
    need not be equally spaced. Returns ``(mid_strikes, densities)`` for the
    interior points (the two endpoints have no central second difference).
    """
    n = len(strikes)
    if n < 3:
        raise ValueError("need at least three strikes")
    if len(calls) != n:
        raise ValueError("strikes and calls must be the same length")
    disc_growth = math.exp(r * t)

    mids, dens = [], []
    for i in range(1, n - 1):
        k0, k1, k2 = strikes[i - 1], strikes[i], strikes[i + 1]
        c0, c1, c2 = calls[i - 1], calls[i], calls[i + 1]
        h1 = k1 - k0
        h2 = k2 - k1
        if h1 <= 0 or h2 <= 0:
            raise ValueError("strikes must be strictly increasing")
        # Non-uniform second derivative (exact for a quadratic through the 3 pts).
        d2 = 2.0 * (h1 * c2 - (h1 + h2) * c1 + h2 * c0) / (h1 * h2 * (h1 + h2))
        pdf = disc_growth * d2
        mids.append(k1)
        dens.append(max(pdf, 0.0))  # clip tiny negative noise
    return mids, dens


def risk_neutral_cdf(strikes: Sequence[float], calls: Sequence[float],
                     t: float, r: float) -> Tuple[List[float], List[float]]:
    """Estimate the risk-neutral CDF F(K) = 1 + e^{rt} dC/dK at midpoints.

    Uses central first differences; returns ``(mid_strikes, cdf_values)``.
    """
    n = len(strikes)
    if n < 3:
        raise ValueError("need at least three strikes")
    disc_growth = math.exp(r * t)
    mids, cdf = [], []
    for i in range(1, n - 1):
        dcdk = (calls[i + 1] - calls[i - 1]) / (strikes[i + 1] - strikes[i - 1])
        val = 1.0 + disc_growth * dcdk
        mids.append(strikes[i])
        cdf.append(min(max(val, 0.0), 1.0))
    return mids, cdf


def price_from_density(strikes: Sequence[float], calls: Sequence[float],
                       t: float, r: float, payoff: Callable[[float], float]) -> float:
    """Price a European payoff by integrating it against the extracted density.

    ``payoff`` maps a terminal underlying value to its cash payoff. The integral
    uses the trapezoidal rule on the recovered density grid; the result is
    discounted at ``r``.
    """
    mids, dens = risk_neutral_density(strikes, calls, t, r)
    if len(mids) < 2:
        raise ValueError("need more strikes to integrate")
    disc = math.exp(-r * t)
    total = 0.0
    for i in range(len(mids) - 1):
        k0, k1 = mids[i], mids[i + 1]
        g0 = payoff(k0) * dens[i]
        g1 = payoff(k1) * dens[i + 1]
        total += 0.5 * (g0 + g1) * (k1 - k0)
    return disc * total


def density_total_mass(strikes: Sequence[float], calls: Sequence[float],
                       t: float, r: float) -> float:
    """Integrate the extracted density; should be close to 1 for a good curve."""
    mids, dens = risk_neutral_density(strikes, calls, t, r)
    total = 0.0
    for i in range(len(mids) - 1):
        total += 0.5 * (dens[i] + dens[i + 1]) * (mids[i + 1] - mids[i])
    return total
