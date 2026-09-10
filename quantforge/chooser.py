"""Simple chooser option (Rubinstein 1991), closed form.

A simple chooser lets the holder decide at a future date ``t_choose`` whether
the option is a call or a put, both struck at ``K`` and expiring at ``T``. At
the choice date the holder picks the more valuable of the two, so by put-call
parity the payoff decomposes into a call struck at ``K`` expiring at ``T`` plus
a put struck at ``K e^{-(r-q)(T - t_choose)}`` expiring at ``t_choose``:

    chooser = C(S, K, T) + P(S, K e^{-(b)(T - t_c)}, t_c)

with cost of carry ``b = r - q``. This is exact; as ``t_choose -> T`` it
approaches a straddle, and as ``t_choose -> 0`` it approaches ``max(call, put)``
today.
"""

import math

from .bsm import (
    call_price, put_price, _validate, delta as bsm_delta, gamma as bsm_gamma,
    vega as bsm_vega, OptionType,
)


def chooser_option(S, K, t_choose, T, r, sigma, b=None) -> float:
    """Price a simple chooser option (Rubinstein 1991).

    Args:
        t_choose: time (years) until the call/put choice is made.
        T: total time (years) to the underlying option's expiry (>= t_choose).
        b: cost of carry (defaults to r).
    """
    _validate(S, K, T, sigma)
    if not (0.0 <= t_choose <= T):
        raise ValueError("require 0 <= t_choose <= T")
    if b is None:
        b = r

    # Long a call to expiry T, plus a put struck at the discounted-forward level
    # expiring at the choice date.
    call = call_price(S, K, T, r, sigma, b=b)
    k_put = K * math.exp(-b * (T - t_choose))
    put = put_price(S, k_put, t_choose, r, sigma, b=b)
    return call + put


def chooser_option_greeks(S, K, t_choose, T, r, sigma, b=None):
    """Greeks of a simple chooser option, exact by decomposition.

    The chooser is exactly ``C(S, K, T) + P(S, K e^{-b(T - t_choose)}, t_choose)``
    -- a call to ``T`` plus a put struck at the discounted-forward level expiring
    at the choice date. Both legs are Black-Scholes prices in ``S`` and ``sigma``
    (the put's strike does not depend on either), so ``delta``, ``gamma``, and
    ``vega`` are the exact sums of the two legs' BSM Greeks -- no finite
    difference. Returns a dict with ``price``, ``delta``, ``gamma``, ``vega``.
    """
    _validate(S, K, T, sigma)
    if not (0.0 <= t_choose <= T):
        raise ValueError("require 0 <= t_choose <= T")
    if b is None:
        b = r
    k_put = K * math.exp(-b * (T - t_choose))

    call = call_price(S, K, T, r, sigma, b=b)
    put = put_price(S, k_put, t_choose, r, sigma, b=b)

    delta = (bsm_delta(S, K, T, r, sigma, OptionType.CALL, b=b)
             + bsm_delta(S, k_put, t_choose, r, sigma, OptionType.PUT, b=b))
    gamma = (bsm_gamma(S, K, T, r, sigma, b=b)
             + bsm_gamma(S, k_put, t_choose, r, sigma, b=b))
    vega = (bsm_vega(S, K, T, r, sigma, b=b)
            + bsm_vega(S, k_put, t_choose, r, sigma, b=b))
    return {"price": call + put, "delta": delta, "gamma": gamma, "vega": vega}
