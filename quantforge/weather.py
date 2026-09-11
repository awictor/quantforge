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
