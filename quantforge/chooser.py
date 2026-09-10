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

from .bsm import call_price, put_price, _validate


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
