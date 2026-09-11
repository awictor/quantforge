"""Weather derivatives: heating/cooling degree days and temperature-index options.

Weather contracts settle on an index built from daily temperatures relative to a
base (conventionally 65 degrees F / 18 degrees C):

    HDD_day = max(base - T_avg, 0)     (heating demand on cold days)
    CDD_day = max(T_avg - base, 0)     (cooling demand on hot days)

A seasonal contract accumulates these over the period and pays a fixed tick value
per index point, often with a cap. This module accumulates degree days, values
the linear (swap) payoff, and prices options on the accumulated index with the
Bachelier (normal) model -- appropriate because a degree-day total is a sum of
many daily contributions and is well approximated as normal. Pure standard
library.
"""

import math


def heating_degree_days(temps, base=65.0):
    """Accumulated heating degree days ``sum_d max(base - T_d, 0)`` over the period.

    ``temps`` is the sequence of daily average temperatures. Each cold day (below
    ``base``) contributes its shortfall; warm days contribute nothing.
    """
    return sum(max(base - t, 0.0) for t in temps)


def cooling_degree_days(temps, base=65.0):
    """Accumulated cooling degree days ``sum_d max(T_d - base, 0)`` over the period."""
    return sum(max(t - base, 0.0) for t in temps)


def degree_day_index(temps, base=65.0, kind="HDD"):
    """Accumulated degree-day index of the requested ``kind`` ("HDD" or "CDD")."""
    if kind == "HDD":
        return heating_degree_days(temps, base)
    if kind == "CDD":
        return cooling_degree_days(temps, base)
    raise ValueError("kind must be 'HDD' or 'CDD'")


def degree_day_swap_payoff(index, strike, tick_value, notional_side=1.0):
    """Linear (swap) payoff on a degree-day index: ``side * tick * (index - strike)``.

    ``tick_value`` is the currency amount per index point; ``notional_side`` is
    ``+1`` for the long-index side (gains when the index exceeds the strike) and
    ``-1`` for the short.
    """
    return notional_side * tick_value * (index - strike)


def degree_day_swap_rate(expected_index):
    """Fair fixed strike of a degree-day swap: the expected accumulated index.

    A degree-day swap pays ``tick * (index - strike)``; its expected value is zero
    when the strike equals the expected index, so the fair strike is
    ``expected_index`` itself.
    """
    return expected_index


def degree_day_collar(expected_index, cap_strike, floor_strike, sigma, r,
                      expiry, tick_value):
    """Zero-cost-style degree-day collar: long a call, short a put.

    Buys protection above ``cap_strike`` (a call) and finances it by selling a put
    struck at ``floor_strike``. Value is
    ``degree_day_option(call, K=cap) - degree_day_option(put, K=floor)``. When both
    strikes coincide the collar reduces to the discounted forward payoff
    ``e^{-r T} tick (expected_index - strike)`` by put-call parity.
    """
    if floor_strike > cap_strike:
        raise ValueError("floor_strike must not exceed cap_strike")
    call = degree_day_option(expected_index, cap_strike, sigma, r, expiry,
                             tick_value, is_call=True)
    put = degree_day_option(expected_index, floor_strike, sigma, r, expiry,
                            tick_value, is_call=False)
    return call - put


def degree_day_option_mc(daily_means, daily_sigma, base, strike, r, expiry,
                         tick_value, kind="HDD", is_call=True, n_paths=20000,
                         seed=4321):
    """Monte Carlo degree-day option over simulated daily temperatures.

    Simulates each day's average temperature as independent normal
    ``N(daily_means[d], daily_sigma^2)``, accumulates the HDD/CDD index over the
    period, and averages the discounted option payoff. An independent reference for
    the Bachelier :func:`degree_day_option` (which approximates the accumulated
    index as normal). Deterministic per seed.
    """
    if daily_sigma < 0:
        raise ValueError("daily_sigma must be non-negative")
    if expiry < 0:
        raise ValueError("expiry must be non-negative")
    if n_paths < 1:
        raise ValueError("n_paths must be positive")
    if kind not in ("HDD", "CDD"):
        raise ValueError("kind must be 'HDD' or 'CDD'")
    disc = math.exp(-r * expiry)
    state = seed & 0xFFFFFFFF
    def _unif():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return (state + 0.5) / 0x80000000
    total = 0.0
    n_days = len(daily_means)
    for _ in range(n_paths):
        index = 0.0
        d = 0
        while d < n_days:
            u1 = _unif()
            u2 = _unif()
            rmag = math.sqrt(-2.0 * math.log(u1))
            z1 = rmag * math.cos(2.0 * math.pi * u2)
            z2 = rmag * math.sin(2.0 * math.pi * u2)
            for z in (z1, z2):
                if d >= n_days:
                    break
                temp = daily_means[d] + daily_sigma * z
                if kind == "HDD":
                    index += max(base - temp, 0.0)
                else:
                    index += max(temp - base, 0.0)
                d += 1
        payoff = max(index - strike, 0.0) if is_call else max(strike - index, 0.0)
        total += payoff
    return disc * tick_value * total / n_paths


def degree_day_option(expected_index, strike, sigma, r, expiry, tick_value,
                      is_call=True, cap=None):
    """Bachelier price of an option on an accumulated degree-day index.

    The seasonal degree-day total is modelled as normal with mean
    ``expected_index`` and standard deviation ``sigma`` (in index points), so a
    call (protection against a high index) or put (low index) is priced by the
    Bachelier formula and scaled by ``tick_value``. With ``m = expected_index -
    strike`` for a call (``strike - expected_index`` for a put) and ``s = sigma``:

        undiscounted = m Phi(m/s) + s phi(m/s)
        price = e^{-r T} * tick_value * undiscounted

    An optional ``cap`` limits the maximum payoff (points), pricing the capped leg
    as a call spread: ``value(strike) - value(strike + cap)``. Put and call satisfy
    ``C - P = e^{-r T} tick (expected_index - strike)`` when uncapped.
    """
    from .mathfns import norm_cdf
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    if expiry < 0:
        raise ValueError("expiry must be non-negative")
    if cap is not None and cap <= 0:
        raise ValueError("cap must be positive")
    disc = math.exp(-r * expiry)

    def _bachelier(k):
        m = (expected_index - k) if is_call else (k - expected_index)
        if sigma == 0.0 or expiry == 0.0:
            return max(m, 0.0)
        d = m / sigma
        pdf = math.exp(-0.5 * d * d) / math.sqrt(2.0 * math.pi)
        return m * norm_cdf(d) + sigma * pdf

    if cap is None:
        undiscounted = _bachelier(strike)
    elif is_call:
        # Capped call = long call(K) - short call(K + cap).
        undiscounted = _bachelier(strike) - _bachelier(strike + cap)
    else:
        # Capped put = long put(K) - short put(K - cap).
        undiscounted = _bachelier(strike) - _bachelier(strike - cap)
    return disc * tick_value * undiscounted
