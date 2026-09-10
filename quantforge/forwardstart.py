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
