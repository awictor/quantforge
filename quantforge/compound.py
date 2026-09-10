"""Compound options (Geske 1979): an option on an option, closed form.

A compound option gives the right, at ``t1``, to pay a first strike ``K1`` for
an underlying option struck at ``K2`` expiring at ``t2 > t1``. The four kinds
are call-on-call, call-on-put, put-on-call, put-on-put.

Geske's formula prices these using the bivariate normal CDF with correlation
``sqrt(t1 / t2)``. It first solves for the critical spot ``S*`` at which the
underlying option is worth exactly ``K1`` at ``t1`` (the exercise boundary of
the compound), then combines univariate and bivariate normal terms.

The bivariate normal CDF is reused from :mod:`quantforge.american`. As ``K1 ->
0`` a call-on-call collapses to the underlying vanilla call, which the tests
check.
"""

import math

from .mathfns import norm_cdf
from .bsm import call_price, put_price, OptionType, _coerce_type, _validate
from .american import _bivariate_normal


def _critical_spot(K1, K2, t2, t1, r, sigma, b, underlying_is_call):
    """Spot S* at t1 where the underlying option is worth exactly K1.

    Solved by bisection; the underlying value is monotincreasing in S for a
    call (decreasing for a put), so the root is unique.
    """
    tau = t2 - t1

    def uval(S):
        if underlying_is_call:
            return call_price(S, K2, tau, r, sigma, b=b)
        return put_price(S, K2, tau, r, sigma, b=b)

    lo, hi = 1e-8, max(K2, 1.0)
    if underlying_is_call:
        # Call value grows with S; expand hi until it exceeds K1.
        while uval(hi) < K1:
            hi *= 2.0
            if hi > 1e12:
                break
    else:
        # Put value falls with S. At hi=K2 the put may still be worth more than
        # K1, so the root S* sits above K2 -> expand hi upward until the put
        # value drops below K1 (bracketing the root between lo and hi).
        while uval(hi) > K1:
            hi *= 2.0
            if hi > 1e12:
                break
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        v = uval(mid) - K1
        if abs(v) < 1e-10:
            return mid
        # Call value increases in S; put decreases.
        if (v > 0) == (not underlying_is_call):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def compound_option(S, K1, K2, t1, t2, r, sigma,
                    kind="call-on-call", b=None) -> float:
    """Price a compound option (Geske 1979).

    Args:
        K1: strike of the compound (paid at ``t1`` to obtain the underlying).
        K2: strike of the underlying option (expiring at ``t2``).
        t1: expiry of the compound (the decision date), ``0 < t1 < t2``.
        t2: expiry of the underlying option.
        kind: "call-on-call", "call-on-put", "put-on-call", "put-on-put".
        b: cost of carry (defaults to r).
    """
    _validate(S, K2, t2, sigma)
    if not (0.0 < t1 < t2):
        raise ValueError("require 0 < t1 < t2")
    if b is None:
        b = r
    kind = kind.lower()
    if kind not in ("call-on-call", "call-on-put", "put-on-call", "put-on-put"):
        raise ValueError("unknown compound kind")

    underlying_is_call = kind.endswith("call")
    S_star = _critical_spot(K1, K2, t2, t1, r, sigma, b, underlying_is_call)

    v1 = sigma * math.sqrt(t1)
    v2 = sigma * math.sqrt(t2)
    carry1 = math.exp((b - r) * t1)
    carry2 = math.exp((b - r) * t2)
    disc1 = math.exp(-r * t1)
    disc2 = math.exp(-r * t2)
    rho = math.sqrt(t1 / t2)

    a1 = (math.log(S / S_star) + (b + 0.5 * sigma * sigma) * t1) / v1
    a2 = a1 - v1
    y1 = (math.log(S / K2) + (b + 0.5 * sigma * sigma) * t2) / v2
    y2 = y1 - v2

    cbnd = _bivariate_normal
    if kind == "call-on-call":
        return (S * carry2 * cbnd(a1, y1, rho)
                - K2 * disc2 * cbnd(a2, y2, rho)
                - K1 * disc1 * norm_cdf(a2))
    if kind == "put-on-call":
        return (K2 * disc2 * cbnd(-a2, y2, -rho)
                - S * carry2 * cbnd(-a1, y1, -rho)
                + K1 * disc1 * norm_cdf(-a2))
    if kind == "call-on-put":
        return (K2 * disc2 * cbnd(-a2, -y2, rho)
                - S * carry2 * cbnd(-a1, -y1, rho)
                - K1 * disc1 * norm_cdf(-a2))
    # put-on-put
    return (S * carry2 * cbnd(a1, -y1, -rho)
            - K2 * disc2 * cbnd(a2, -y2, -rho)
            + K1 * disc1 * norm_cdf(a2))
