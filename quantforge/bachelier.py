"""Bachelier (normal) model: pricing, Greeks, and implied normal volatility.

The Bachelier model assumes the underlying (or forward) follows *arithmetic*
Brownian motion rather than geometric:

    dF = sigma_N dW

so the terminal value is normally distributed and can go negative. This is the
market standard for interest-rate and spread options, where forwards routinely
trade at or below zero and a lognormal model breaks down. ``sigma_N`` is the
**normal** (absolute, price-unit) volatility, not a percentage vol.

The undiscounted call value on a forward ``F`` struck at ``K`` is

    C = e^{-r t} [ (F - K) N(d) + sigma_N sqrt(t) n(d) ],   d = (F - K)/(sigma_N sqrt(t))

with the put following from parity. All functions take the forward ``F``
directly; for a spot price ``S`` with carry ``b`` use ``F = S e^{b t}``.
"""

import math

from .mathfns import norm_cdf, norm_pdf
from .bsm import OptionType, _coerce_type


def _validate(F, K, t, sigma):
    if t < 0:
        raise ValueError("t must be non-negative")
    if sigma < 0:
        raise ValueError("normal volatility must be non-negative")


def bachelier_price(F, K, t, r, sigma, option_type=OptionType.CALL) -> float:
    """Bachelier price of a European option on a forward ``F``.

    ``sigma`` is the normal (absolute) volatility. ``r`` discounts the payoff
    from expiry; pass ``r=0`` to price on the forward directly.
    """
    ot = _coerce_type(option_type)
    _validate(F, K, t, sigma)
    disc = math.exp(-r * t)

    if t == 0 or sigma == 0:
        payoff = max(F - K, 0.0) if ot is OptionType.CALL else max(K - F, 0.0)
        return disc * payoff

    vsqrt = sigma * math.sqrt(t)
    d = (F - K) / vsqrt
    if ot is OptionType.CALL:
        return disc * ((F - K) * norm_cdf(d) + vsqrt * norm_pdf(d))
    return disc * ((K - F) * norm_cdf(-d) + vsqrt * norm_pdf(d))


def bachelier_delta(F, K, t, r, sigma, option_type=OptionType.CALL) -> float:
    """dPrice/dF (in the forward). Call delta is e^{-rt} N(d)."""
    ot = _coerce_type(option_type)
    _validate(F, K, t, sigma)
    disc = math.exp(-r * t)
    if t == 0 or sigma == 0:
        itm = F > K if ot is OptionType.CALL else F < K
        base = disc if itm else 0.0
        return base if ot is OptionType.CALL else -base
    d = (F - K) / (sigma * math.sqrt(t))
    return disc * norm_cdf(d) if ot is OptionType.CALL else -disc * norm_cdf(-d)


def bachelier_gamma(F, K, t, r, sigma) -> float:
    """d2Price/dF2. Same for calls and puts."""
    _validate(F, K, t, sigma)
    if t == 0 or sigma == 0:
        return 0.0
    vsqrt = sigma * math.sqrt(t)
    d = (F - K) / vsqrt
    return math.exp(-r * t) * norm_pdf(d) / vsqrt


def bachelier_vega(F, K, t, r, sigma) -> float:
    """dPrice/dsigma_N (per unit of normal vol). Same for calls and puts."""
    _validate(F, K, t, sigma)
    if t == 0 or sigma == 0:
        return 0.0
    vsqrt = sigma * math.sqrt(t)
    d = (F - K) / vsqrt
    return math.exp(-r * t) * math.sqrt(t) * norm_pdf(d)


def bachelier_theta(F, K, t, r, sigma, option_type=OptionType.CALL) -> float:
    """Calendar theta ``-dPrice/dt`` in the Bachelier (normal) model, analytic.

    Differentiating the discounted normal price gives, for a call,

        theta = r * price - e^{-rt} [ sigma phi(d) / (2 sqrt(t)) ],

    where ``d = (F - K)/(sigma sqrt(t))``. The normal-model time value grows
    with maturity (the ``sigma phi/(2 sqrt t)`` piece, same for calls and puts,
    enters ``-dP/dt`` with a minus sign); the ``r * price`` piece is the
    discount drift. At ``r = 0`` this is the pure decay
    ``-sigma phi(d)/(2 sqrt(t))``.
    """
    ot = _coerce_type(option_type)
    _validate(F, K, t, sigma)
    if t == 0 or sigma == 0:
        return 0.0
    disc = math.exp(-r * t)
    vsqrt = sigma * math.sqrt(t)
    d = (F - K) / vsqrt
    price = bachelier_price(F, K, t, r, sigma, ot)
    decay = disc * sigma * norm_pdf(d) / (2.0 * math.sqrt(t))
    return r * price - decay


def bachelier_cash_or_nothing(F, K, t, r, sigma, option_type=OptionType.CALL,
                              cash=1.0) -> float:
    """Bachelier cash-or-nothing digital: pays ``cash`` if in the money.

    In the normal model ``F_T`` is Gaussian, so with ``d = (F - K)/(sigma sqrt t)``
    a call (pays when ``F_T > K``) is ``cash e^{-rt} N(d)`` and a put is
    ``cash e^{-rt} N(-d)``.
    """
    ot = _coerce_type(option_type)
    _validate(F, K, t, sigma)
    disc = math.exp(-r * t)
    if t == 0 or sigma == 0:
        itm = F > K if ot is OptionType.CALL else F < K
        return disc * cash if itm else 0.0
    d = (F - K) / (sigma * math.sqrt(t))
    return cash * disc * (norm_cdf(d) if ot is OptionType.CALL else norm_cdf(-d))


def bachelier_asset_or_nothing(F, K, t, r, sigma,
                               option_type=OptionType.CALL) -> float:
    """Bachelier asset-or-nothing digital: pays the forward ``F_T`` if in the money.

    With ``F_T`` Gaussian, ``E[F_T 1_{F_T > K}] = F N(d) + sigma sqrt(t) phi(d)``
    (call), discounted at ``r``. Note the vanilla Bachelier call equals this
    asset-or-nothing minus ``K`` times the cash-or-nothing, mirroring the
    Black-Scholes decomposition.
    """
    ot = _coerce_type(option_type)
    _validate(F, K, t, sigma)
    disc = math.exp(-r * t)
    if t == 0 or sigma == 0:
        itm = F > K if ot is OptionType.CALL else F < K
        return disc * F if itm else 0.0
    vsqrt = sigma * math.sqrt(t)
    d = (F - K) / vsqrt
    if ot is OptionType.CALL:
        return disc * (F * norm_cdf(d) + vsqrt * norm_pdf(d))
    return disc * (F * norm_cdf(-d) - vsqrt * norm_pdf(d))


def bachelier_implied_vol(target_price, F, K, t, r, option_type=OptionType.CALL,
                          tol=1e-10, max_iter=100) -> float:
    """Solve for the normal volatility that reproduces ``target_price``.

    Newton's method on vega with a bisection fallback. Rejects prices outside
    the no-arbitrage band [intrinsic, forward-bound].
    """
    ot = _coerce_type(option_type)
    if t <= 0:
        raise ValueError("cannot imply vol at or past expiry")
    disc = math.exp(-r * t)
    intrinsic = disc * (max(F - K, 0.0) if ot is OptionType.CALL else max(K - F, 0.0))
    if target_price < intrinsic - 1e-12:
        raise ValueError("price below intrinsic value")

    # Upper vol bound: grow until the model price exceeds the target.
    hi = max(abs(F - K), 1.0)
    while bachelier_price(F, K, t, r, hi, ot) < target_price:
        hi *= 2.0
        if hi > 1e12:
            raise ValueError("price above the no-arbitrage upper bound")
    lo = 0.0

    # Seed from the ATM inversion (exact at F==K): sigma ~ price*sqrt(2pi/t).
    # Off-ATM this is only a starting guess; the maintained [lo, hi] bracket
    # guarantees convergence regardless.
    sigma = (target_price / disc) * math.sqrt(2.0 * math.pi / t)
    sigma = min(max(sigma, 1e-8), hi)

    for _ in range(max_iter):
        price = bachelier_price(F, K, t, r, sigma, ot)
        diff = price - target_price
        if abs(diff) < tol:
            return sigma
        if diff > 0:
            hi = sigma
        else:
            lo = sigma
        vega = bachelier_vega(F, K, t, r, sigma)
        if vega > 1e-14:
            step = sigma - diff / vega
            if lo < step < hi:
                sigma = step
                continue
        sigma = 0.5 * (lo + hi)
    return sigma


def bachelier_greeks(F, K, t, r, sigma, option_type=OptionType.CALL):
    """Bundle the Bachelier Greeks: delta, gamma, vega (analytic) plus theta.

    ``delta``, ``gamma``, and ``vega`` reuse the exact closed forms
    :func:`bachelier_delta`, :func:`bachelier_gamma`, :func:`bachelier_vega`;
    ``theta`` (calendar decay, ``-dV/dt``) is a central finite difference of
    :func:`bachelier_price`. Returns a dict with ``price``, ``delta``, ``gamma``,
    ``vega``, ``theta``. All are in normal-model (absolute-vol) terms.
    """
    ot = _coerce_type(option_type)
    _validate(F, K, t, sigma)
    price = bachelier_price(F, K, t, r, sigma, ot)
    delta = bachelier_delta(F, K, t, r, sigma, ot)
    gamma = bachelier_gamma(F, K, t, r, sigma)
    vega = bachelier_vega(F, K, t, r, sigma)
    theta = bachelier_theta(F, K, t, r, sigma, ot)
    return {"price": price, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}
