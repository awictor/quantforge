"""Discrete-monitoring geometric-average Asian (exotics.discrete_geometric_asian)."""

import math

import pytest

from quantforge import (
    discrete_geometric_asian, discrete_geometric_asian_greeks, geometric_asian,
)
from quantforge.bsm import price, OptionType


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2


def test_single_fixing_recovers_vanilla():
    for ot in (OptionType.CALL, OptionType.PUT):
        assert discrete_geometric_asian(S, K, T, R, SIG, n_fixings=1,
                                        option_type=ot) == pytest.approx(
            price(S, K, T, R, SIG, ot), abs=1e-10)


def test_converges_to_continuous():
    cont = geometric_asian(S, K, T, R, SIG, "call")
    p50 = discrete_geometric_asian(S, K, T, R, SIG, n_fixings=50, option_type="call")
    p500 = discrete_geometric_asian(S, K, T, R, SIG, n_fixings=500, option_type="call")
    # Monotone approach from above; the finer grid is strictly closer.
    assert abs(p500 - cont) < abs(p50 - cont)
    assert p500 == pytest.approx(cont, abs=1e-2)


def test_explicit_fixing_times_equal_equally_spaced():
    n = 4
    times = [T * i / n for i in range(1, n + 1)]
    a = discrete_geometric_asian(S, K, T, R, SIG, n_fixings=n, option_type="call")
    b = discrete_geometric_asian(S, K, T, R, SIG, fixing_times=times,
                                 option_type="call")
    assert a == pytest.approx(b, abs=1e-12)


def _mc(n, ot, npaths=300000, seed=3):
    import random
    rng = random.Random(seed)
    acc = 0.0
    for _ in range(npaths):
        logS = math.log(S)
        logsum = 0.0
        prevt = 0.0
        for i in range(1, n + 1):
            ti = T * i / n
            dt = ti - prevt
            prevt = ti
            z = rng.gauss(0.0, 1.0)
            logS += (R - 0.5 * SIG * SIG) * dt + SIG * math.sqrt(dt) * z
            logsum += logS
        G = math.exp(logsum / n)
        pay = max(G - K, 0.0) if ot == "call" else max(K - G, 0.0)
        acc += pay
    return math.exp(-R * T) * acc / npaths


@pytest.mark.parametrize("ot", ["call", "put"])
def test_matches_monte_carlo(ot):
    n = 12
    cf = discrete_geometric_asian(S, K, T, R, SIG, n_fixings=n, option_type=ot)
    mc = _mc(n, ot)
    assert cf == pytest.approx(mc, rel=0.01)


def test_geometric_below_or_equal_arithmetic_intuition():
    # Averaging lowers effective vol, so a discrete geo Asian call is cheaper
    # than the vanilla (single-fixing) call of the same params.
    van = price(S, K, T, R, SIG, "call")
    asian = discrete_geometric_asian(S, K, T, R, SIG, n_fixings=12, option_type="call")
    assert asian < van


def test_greeks_signs_and_price_field():
    g = discrete_geometric_asian_greeks(S, K, T, R, SIG, n_fixings=12,
                                        option_type="call")
    assert g["delta"] > 0.0
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0
    assert g["price"] == pytest.approx(
        discrete_geometric_asian(S, K, T, R, SIG, n_fixings=12, option_type="call"),
        abs=1e-9)


def test_validation():
    with pytest.raises(ValueError):
        discrete_geometric_asian(S, K, T, R, SIG)  # neither n nor times
    with pytest.raises(ValueError):
        discrete_geometric_asian(S, K, T, R, SIG, fixing_times=[0.5, 1.5])  # > t
