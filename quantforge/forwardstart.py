"""Forward-start and cliquet option pricing (closed form under BSM).

A **forward-start** option is set live today but its strike is fixed at a future
date ``t_start`` as a multiple ``alpha`` of the then-prevailing spot
(``K = alpha * S_{t_start}``). Under Black-Scholes-Merton with constant
parameters, the price scales with the current spot and does not depend on the
strike level, because the payoff is homogeneous in the (unknown) future spot.
Rubinstein's (1990) result gives the value as the current spot times a scaled
Black-Scholes price on a unit underlying:

    FS = S * e^{(b - r) t_start} * BS(1, alpha, tau, r, sigma, b)

where ``tau = t_expiry - t_start`` is the life of the option after it starts.

A **cliquet** (ratchet) is a strip of consecutive forward-start options that
reset periodically; its value is the sum of the constituent forward-starts.
"""

import math
from typing import Sequence

from .bsm import OptionType, _coerce_type, _validate, price as bsm_price


def forward_start_price(S, t_start, t_expiry, r, sigma, alpha=1.0,
                        option_type=OptionType.CALL, b=None) -> float:
    """Price a forward-start option whose strike is set at ``t_start``.

    Args:
        S: current spot.
        t_start: time (years) until the strike is fixed. 0 reduces to a vanilla.
        t_expiry: total time (years) to the option's expiry (> t_start).
        alpha: moneyness multiple; strike = alpha * S_{t_start}. alpha=1 is ATM.
        b: cost of carry (defaults to r).

    Returns the present value.
    """
    ot = _coerce_type(option_type)
    _validate(S, 1.0, t_expiry, sigma)  # K placeholder; validated separately
    if alpha <= 0:
        raise ValueError("alpha (moneyness multiple) must be positive")
    if t_start < 0 or t_expiry <= 0:
        raise ValueError("need t_start >= 0 and t_expiry > 0")
    if t_start >= t_expiry:
        raise ValueError("t_start must be before t_expiry")
    if b is None:
        b = r

    tau = t_expiry - t_start
    # Value per unit of underlying-at-reset, discounted for the carry earned
    # between now and the reset date.
    unit = bsm_price(1.0, alpha, tau, r, sigma, ot, b=b)
    return S * math.exp((b - r) * t_start) * unit


def forward_start_greeks(S, t_start, t_expiry, r, sigma, alpha=1.0,
                         option_type=OptionType.CALL, b=None):
    """Greeks of a forward-start option (Rubinstein), exact where possible.

    The price is ``FS = S e^{(b-r) t_start} * u`` where ``u`` is a unit
    Black-Scholes price on a unit underlying and does **not** depend on ``S``.
    So the value is exactly linear in the spot: ``delta = e^{(b-r) t_start} u``
    (constant in ``S``) and ``gamma = 0`` -- a forward-start has no spot gamma
    until its strike is fixed. ``vega`` and ``theta`` (calendar decay, both
    ``t_start`` and ``t_expiry`` shifting together) are central finite
    differences of the closed form. Returns a dict with ``price``, ``delta``,
    ``gamma``, ``vega``, ``theta``.
    """
    ot = _coerce_type(option_type)
    if alpha <= 0:
        raise ValueError("alpha (moneyness multiple) must be positive")
    if t_start < 0 or t_expiry <= 0:
        raise ValueError("need t_start >= 0 and t_expiry > 0")
    if t_start >= t_expiry:
        raise ValueError("t_start must be before t_expiry")
    if b is None:
        b = r

    tau = t_expiry - t_start
    unit = bsm_price(1.0, alpha, tau, r, sigma, ot, b=b)
    scale = math.exp((b - r) * t_start)
    price = S * scale * unit
    delta = scale * unit          # exact: price is linear in S
    gamma = 0.0                   # no spot convexity before the strike is set

    def px(sigma_=sigma, dt=0.0):
        return forward_start_price(S, t_start - dt, t_expiry - dt, r, sigma_,
                                   alpha, ot, b=b)

    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t_start) if t_start > 0 else 1e-4
    theta = -(px(dt=-ht) - px(dt=ht)) / (2.0 * ht)
    return {"price": price, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


def cliquet_price(S, reset_times: Sequence[float], r, sigma, alpha=1.0,
                  option_type=OptionType.CALL, b=None) -> float:
    """Price a cliquet (ratchet) as a strip of forward-start options.

    Args:
        reset_times: increasing schedule of reset/expiry dates in years, e.g.
            [0.25, 0.5, 0.75, 1.0]. Each consecutive pair (t_i, t_{i+1}) is one
            forward-start period that starts at t_i and expires at t_{i+1}. The
            first period starts now (t=0) and ends at reset_times[0].

    Returns the total present value of the strip.
    """
    ot = _coerce_type(option_type)
    if not reset_times:
        raise ValueError("cliquet needs at least one reset time")
    times = list(reset_times)
    if any(times[i] >= times[i + 1] for i in range(len(times) - 1)):
        raise ValueError("reset_times must be strictly increasing")
    if times[0] <= 0:
        raise ValueError("first reset time must be positive")
    if b is None:
        b = r

    total = 0.0
    starts = [0.0] + times[:-1]
    for t_start, t_end in zip(starts, times):
        if t_start == 0.0:
            # First period is a plain (spot-strike) option struck at alpha*S.
            total += bsm_price(S, alpha * S, t_end, r, sigma, ot, b=b)
        else:
            total += forward_start_price(S, t_start, t_end, r, sigma, alpha, ot, b=b)
    return total


def cliquet_greeks(S, reset_times: Sequence[float], r, sigma, alpha=1.0,
                   option_type=OptionType.CALL, b=None):
    """Greeks of a cliquet (ratchet) by central finite differences.

    A cliquet is a strip of consecutive forward-start options. Only the first
    (spot-strike) period carries spot gamma; every later forward-start period is
    linear in the current spot, so the cliquet's ``gamma`` comes entirely from
    the first period and is small relative to a single vanilla. ``delta``,
    ``gamma``, ``vega``, and ``theta`` are central finite differences of
    :func:`cliquet_price`; ``theta`` shifts every reset date together. Returns a
    dict with ``price``, ``delta``, ``gamma``, ``vega``, ``theta``.
    """
    ot = _coerce_type(option_type)
    if not reset_times:
        raise ValueError("cliquet needs at least one reset time")
    times = list(reset_times)
    if any(times[i] >= times[i + 1] for i in range(len(times) - 1)):
        raise ValueError("reset_times must be strictly increasing")
    if times[0] <= 0:
        raise ValueError("first reset time must be positive")
    if b is None:
        b = r

    def px(S_=S, sigma_=sigma, dt=0.0):
        shifted = [tm - dt for tm in times]
        return cliquet_price(S_, shifted, r, sigma_, alpha, ot, b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * times[0])
    theta = -(px(dt=-ht) - px(dt=ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}
