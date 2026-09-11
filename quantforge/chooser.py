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
from .american import _bivariate_normal
from .mathfns import norm_cdf


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


def _chooser_critical_spot(Kc, Kp, tc, Tc, Tp, r, sigma, b):
    """Solve for the spot ``I`` at the choice date where the call and put legs
    are equally valuable: ``C(I, Kc, Tc - tc) = P(I, Kp, Tp - tc)``.

    The call value rises and the put value falls in the spot, so their
    difference is monotone increasing and a bisection has a unique root.
    """
    tau_c = Tc - tc
    tau_p = Tp - tc

    def diff(x):
        return (call_price(x, Kc, tau_c, r, sigma, b=b)
                - put_price(x, Kp, tau_p, r, sigma, b=b))

    lo, hi = 1e-8, max(Kc, Kp)
    # Expand the upper bracket until the difference turns positive.
    while diff(hi) < 0.0:
        hi *= 2.0
        if hi > 1e12:
            break
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if diff(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def complex_chooser_option(S, Kc, Kp, t_choose, Tc, Tp, r, sigma, b=None) -> float:
    """Complex chooser option (Rubinstein 1991), closed form.

    At the choice date ``t_choose`` the holder keeps whichever is worth more: a
    call struck at ``Kc`` expiring at ``Tc``, or a put struck at ``Kp`` expiring
    at ``Tp`` (the two legs may differ in both strike and maturity). Rubinstein's
    formula prices this with bivariate normals coupling the choice date to each
    leg's expiry:

        V = S e^{(b-r)Tc} M(d1, y1; rho_c) - Kc e^{-r Tc} M(d2, y1 - sig sqrt Tc; rho_c)
            - S e^{(b-r)Tp} M(-d1, -y2; rho_p) + Kp e^{-r Tp} M(-d2, -y2 + sig sqrt Tp; rho_p)

    where ``d1,d2`` use the critical spot ``I`` (the level where the two legs are
    equal at ``t_choose``), ``y1,y2`` use ``Kc,Kp``, and
    ``rho_c = sqrt(t_choose/Tc)``, ``rho_p = sqrt(t_choose/Tp)``.
    """
    if b is None:
        b = r
    _validate(S, Kc, Tc, sigma)
    _validate(S, Kp, Tp, sigma)
    if not (0.0 < t_choose <= min(Tc, Tp)):
        raise ValueError("require 0 < t_choose <= min(Tc, Tp)")

    I = _chooser_critical_spot(Kc, Kp, t_choose, Tc, Tp, r, sigma, b)
    st = sigma * math.sqrt(t_choose)
    d1 = (math.log(S / I) + (b + 0.5 * sigma * sigma) * t_choose) / st
    d2 = d1 - st
    y1 = (math.log(S / Kc) + (b + 0.5 * sigma * sigma) * Tc) / (sigma * math.sqrt(Tc))
    y2 = (math.log(S / Kp) + (b + 0.5 * sigma * sigma) * Tp) / (sigma * math.sqrt(Tp))
    rho_c = math.sqrt(t_choose / Tc)
    rho_p = math.sqrt(t_choose / Tp)
    cc = math.exp((b - r) * Tc)
    cp = math.exp((b - r) * Tp)
    return (
        S * cc * _bivariate_normal(d1, y1, rho_c)
        - Kc * math.exp(-r * Tc) * _bivariate_normal(d2, y1 - sigma * math.sqrt(Tc), rho_c)
        - S * cp * _bivariate_normal(-d1, -y2, rho_p)
        + Kp * math.exp(-r * Tp) * _bivariate_normal(-d2, -y2 + sigma * math.sqrt(Tp), rho_p)
    )


def complex_chooser_option_greeks(S, Kc, Kp, t_choose, Tc, Tp, r, sigma, b=None):
    """Greeks of a complex chooser by central finite differences of
    :func:`complex_chooser_option`: ``delta`` (dV/dS), ``gamma`` (d2V/dS2),
    ``vega`` (dV/dsigma), ``theta`` (calendar decay -- all maturities and the
    choice date shrink together). Returns a dict with ``price`` and those
    fields.
    """
    if b is None:
        b = r

    def px(S_=S, sigma_=sigma, shift=0.0):
        return complex_chooser_option(S_, Kc, Kp, t_choose - shift,
                                      Tc - shift, Tp - shift, r, sigma_, b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t_choose)
    theta = -(px(shift=ht) - px(shift=-ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}
