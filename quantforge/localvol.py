"""Dupire local volatility from a call-price surface.

Dupire (1994) showed that a single deterministic local-volatility function
``sigma_loc(K, T)`` reproduces every European call price on the surface. In
terms of the (undiscounted-payoff) call price ``C(K, T)`` it is

    sigma_loc^2(K, T) =
        [ dC/dT + (r - q) K dC/dK + q C ] / ( 0.5 K^2 d2C/dK2 ).

Given a callable ``C(K, T)`` (from a fitted implied-vol surface, say) we
evaluate the three derivatives by central finite differences. Under
Black-Scholes with a flat implied vol the local vol equals that constant,
which the tests use as a check.
"""

import math
from typing import Callable


def dupire_local_vol(call_fn: Callable[[float, float], float],
                     K: float, T: float, r: float, q: float = 0.0,
                     dK: float = None, dT: float = None) -> float:
    """Dupire local volatility at strike ``K`` and maturity ``T``.

    Args:
        call_fn: ``C(K, T)`` returning the European call price for strike K and
            maturity T (both positive). Must be evaluable in a neighborhood of
            (K, T) for the finite differences.
        r, q: risk-free rate and continuous dividend yield.
        dK, dT: finite-difference bumps; default to small fractions of K and T.

    Returns the local volatility (not variance). Raises if the local variance
    comes out non-positive (a sign of an arbitrageable / too-noisy surface).
    """
    if K <= 0 or T <= 0:
        raise ValueError("K and T must be positive")
    if dK is None:
        dK = max(1e-3, 1e-3 * K)
    if dT is None:
        dT = max(1e-4, 1e-3 * T)
    # Keep the maturity bump inside (0, T].
    dT = min(dT, 0.5 * T)

    C = call_fn(K, T)
    # First derivative in maturity (central, or forward near T small).
    C_Tup = call_fn(K, T + dT)
    C_Tdn = call_fn(K, T - dT)
    dC_dT = (C_Tup - C_Tdn) / (2.0 * dT)

    # Strike derivatives (central).
    C_Kup = call_fn(K + dK, T)
    C_Kdn = call_fn(K - dK, T)
    dC_dK = (C_Kup - C_Kdn) / (2.0 * dK)
    d2C_dK2 = (C_Kup - 2.0 * C + C_Kdn) / (dK * dK)

    numerator = dC_dT + (r - q) * K * dC_dK + q * C
    denominator = 0.5 * K * K * d2C_dK2

    if denominator <= 0:
        raise ValueError(
            "non-positive convexity (d2C/dK2 <= 0); surface is arbitrageable "
            "or the strike bump is too small/noisy"
        )
    local_var = numerator / denominator
    if local_var <= 0:
        raise ValueError(
            f"non-positive local variance ({local_var:.3g}); check the surface "
            "or the finite-difference bumps"
        )
    return math.sqrt(local_var)


def local_vol_from_implied(implied_vol_fn: Callable[[float, float], float],
                           S: float, K: float, T: float, r: float, q: float = 0.0,
                           dK: float = None, dT: float = None) -> float:
    """Dupire local vol from an implied-vol surface ``sigma_imp(K, T)``.

    Wraps :func:`dupire_local_vol` by turning the implied-vol surface into a
    call-price surface with the Black-Scholes-Merton formula (carry ``b = r-q``).
    """
    from .bsm import call_price

    def call_fn(k, t):
        sig = implied_vol_fn(k, t)
        return call_price(S, k, t, r, sig, b=r - q)

    return dupire_local_vol(call_fn, K, T, r, q, dK, dT)
