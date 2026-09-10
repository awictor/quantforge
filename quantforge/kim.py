"""American option pricing via Kim's (1990) integral equation.

An American option's value equals the matching European value plus an
early-exercise premium, an integral over the (unknown) early-exercise boundary
``B(t)``. Kim (1990) writes that boundary as the solution of a Volterra integral
equation: at each time-to-maturity the value-matching condition (the option is
worth its intrinsic value on the boundary) becomes

    K - B(tau) = P_E(B(tau), tau)
        + integral_0^tau [ r K e^{-r s} N(-d2) - q B(tau) e^{-q s} N(-d1) ] ds

for an American put, with ``d1, d2`` built from ``B(tau)`` and ``B(tau - s)``.
We discretise time, march from expiry (where the boundary is known in closed
form) backward, and solve the one-dimensional value-matching equation for the
boundary at each step by bisection. The price is then the European value plus
the premium integral evaluated at the current spot. Pure standard library; the
converged boundary and price agree with a fine binomial tree.
"""

import math

from .bsm import OptionType, _coerce_type, put_price, call_price
from .mathfns import norm_cdf


def _d12(S, B, r, q, sigma, s):
    """Black-Scholes d1, d2 for spot ``S`` against level ``B`` over horizon ``s``."""
    vs = sigma * math.sqrt(s)
    d1 = (math.log(S / B) + (r - q + 0.5 * sigma * sigma) * s) / vs
    return d1, d1 - vs


def _integrand(S, Bs, K, r, q, sigma, s):
    """Kim premium integrand at horizon ``s`` (with its finite s -> 0 limit).

    As ``s -> 0`` with ``S`` strictly above ``Bs`` the integrand vanishes, but at
    ``S == Bs`` (the boundary equation) ``d1, d2 -> 0`` so the limit is
    ``r K / 2 - q S / 2`` (from ``N(0) = 1/2``), NOT zero -- getting this endpoint
    wrong biases the whole premium.
    """
    if s <= 0.0:
        if abs(math.log(S / Bs)) < 1e-12:
            return 0.5 * r * K - 0.5 * q * S
        return 0.0
    d1, d2 = _d12(S, Bs, r, q, sigma, s)
    return (r * K * math.exp(-r * s) * norm_cdf(-d2)
            - q * S * math.exp(-q * s) * norm_cdf(-d1))


def _boundary_put(K, t, r, q, sigma, n_steps):
    """Solve the American-put early-exercise boundary B[i] at times i*dt.

    Backward Volterra recursion: B[N] at expiry is K*min(1, r/q) (K if q=0),
    and each earlier B[i] solves the value-matching equation given the already
    -computed later boundary points (trapezoidal premium integral).
    """
    dt = t / n_steps
    # Terminal boundary: for a put, B(T) = K min(1, r/q); q -> 0 gives K.
    B = [0.0] * (n_steps + 1)
    B[n_steps] = K if q <= 1e-14 else K * min(1.0, r / q)

    for i in range(n_steps - 1, -1, -1):
        tau = (n_steps - i) * dt   # time-to-maturity at grid node i

        def value_match(b):
            # LHS: intrinsic minus European put at (b, tau).
            lhs = (K - b) - put_price(b, K, tau, r, sigma, b=r - q)
            # RHS: premium integral over [0, tau] via the trapezoid rule on the
            # already-known boundary B[i+1..n_steps] (B(tau - s) = B at node i+m).
            integ = 0.0
            prev = None
            for m in range(0, n_steps - i + 1):
                s = m * dt
                # At s = 0 the boundary point B[i] itself is the unknown b, so
                # use b (not B[i]) as the level, giving the correct r K/2 limit.
                Bs = b if m == 0 else B[i + m]
                val = _integrand(b, Bs, K, r, q, sigma, s)
                if prev is not None:
                    integ += 0.5 * (prev + val) * dt
                prev = val
            return lhs - integ

        # Bisect for B[i] in (0, K): value_match is monotone in b here.
        lo, hi = 1e-8 * K, K
        flo, fhi = value_match(lo), value_match(hi)
        if flo * fhi > 0:
            B[i] = B[i + 1]   # fall back to continuity if no sign change
            continue
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            fm = value_match(mid)
            if flo * fm <= 0:
                hi, fhi = mid, fm
            else:
                lo, flo = mid, fm
        B[i] = 0.5 * (lo + hi)
    return B, dt


def _premium_put(S, K, t, r, q, sigma, B, dt):
    """Early-exercise premium of an American put at spot ``S`` (trapezoid)."""
    n_steps = len(B) - 1
    integ = 0.0
    prev = None
    for m in range(0, n_steps + 1):
        s = m * dt
        # Integrand needs the boundary s years in the future, i.e. at
        # time-to-maturity (tau - s). With tau = t (full grid), that is B[m]
        # (ttm decreases from tau at m=0 to 0 at m=n_steps).
        Bs = B[m]
        # S is the (interior) spot, strictly above B(0), so the s -> 0 term is 0.
        val = _integrand(S, Bs, K, r, q, sigma, s)
        if prev is not None:
            integ += 0.5 * (prev + val) * dt
        prev = val
    return integ


def kim_american_put(S, K, t, r, sigma, q=0.0, n_steps=80):
    """American put price via Kim's integral equation.

    Solves the early-exercise boundary on an ``n_steps`` time grid, then returns
    the European put plus the early-exercise premium integrated at spot ``S``.
    If ``S`` is at or below the current boundary the option is exercised, so the
    intrinsic value is returned.
    """
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t < 0 or sigma <= 0:
        raise ValueError("need t >= 0 and sigma > 0")
    if t == 0:
        return max(K - S, 0.0)

    B, dt = _boundary_put(K, t, r, q, sigma, n_steps)
    if S <= B[0]:
        return K - S
    euro = put_price(S, K, t, r, sigma, b=r - q)
    prem = _premium_put(S, K, t, r, q, sigma, B, dt)
    return max(euro + prem, K - S)


def kim_american_call(S, K, t, r, sigma, q=0.0, n_steps=80):
    """American call price via the put-call symmetry for American options.

    A dividend-paying American call maps to an American put by the McDonald-
    Schroder symmetry ``C(S, K, r, q) = P(K, S, q, r)`` (spot<->strike,
    rate<->dividend). With ``q = 0`` the call is never exercised early and this
    returns the European call.
    """
    if q <= 1e-14:
        return call_price(S, K, t, r, sigma, b=r)   # no early exercise
    return kim_american_put(K, S, t, q, sigma, q=r, n_steps=n_steps)


def kim_exercise_boundary(K, t, r, sigma, q=0.0, n_steps=80):
    """Return the American-put early-exercise boundary ``B(t_i)`` on the grid.

    ``B[i]`` is the critical spot at time ``i * (t / n_steps)`` below which
    immediate exercise is optimal; ``B[n_steps]`` is the expiry value.
    """
    B, _dt = _boundary_put(K, t, r, q, sigma, n_steps)
    return B
