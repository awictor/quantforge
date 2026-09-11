"""Discrete-monitoring arithmetic-average Asian, Levy match (discrete_arithmetic_asian)."""

import math

import pytest

from quantforge import (
    discrete_arithmetic_asian, discrete_arithmetic_asian_greeks,
    discrete_geometric_asian,
)
from quantforge.bsm import price, OptionType


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2


def test_single_fixing_recovers_vanilla():
    for ot in (OptionType.CALL, OptionType.PUT):
        assert discrete_arithmetic_asian(S, K, T, R, SIG, n_fixings=1,
                                         option_type=ot) == pytest.approx(
            price(S, K, T, R, SIG, ot), abs=1e-10)


def test_arithmetic_above_geometric():
    # AM-GM: arithmetic average >= geometric, so the arithmetic Asian call
    # is worth more than the geometric one at the same params.
    ar = discrete_arithmetic_asian(S, K, T, R, SIG, n_fixings=12, option_type="call")
    ge = discrete_geometric_asian(S, K, T, R, SIG, n_fixings=12, option_type="call")
    assert ar > ge


def test_explicit_fixing_times_equal_equally_spaced():
    n = 4
    times = [T * i / n for i in range(1, n + 1)]
    a = discrete_arithmetic_asian(S, K, T, R, SIG, n_fixings=n, option_type="call")
    b = discrete_arithmetic_asian(S, K, T, R, SIG, fixing_times=times,
                                  option_type="call")
    assert a == pytest.approx(b, abs=1e-12)


def _mc_control_variate(n, ot, npaths=200000, seed=5):
    """Arithmetic-Asian MC using the exact discrete geometric Asian as control."""
    import random
    rng = random.Random(seed)
    ge_exact = discrete_geometric_asian(S, K, T, R, SIG, n_fixings=n, option_type=ot)
    disc = math.exp(-R * T)
    sa = sg = 0.0
    for _ in range(npaths):
        logS = math.log(S)
        prevt = 0.0
        asum = 0.0
        gsum = 0.0
        for i in range(1, n + 1):
            ti = T * i / n
            dt = ti - prevt
            prevt = ti
            z = rng.gauss(0.0, 1.0)
            logS += (R - 0.5 * SIG * SIG) * dt + SIG * math.sqrt(dt) * z
            asum += math.exp(logS)
            gsum += logS
        A = asum / n
        G = math.exp(gsum / n)
        pa = max(A - K, 0.0) if ot == "call" else max(K - A, 0.0)
        pg = max(G - K, 0.0) if ot == "call" else max(K - G, 0.0)
        sa += pa
        sg += pg
    mc_a = disc * sa / npaths
    mc_g = disc * sg / npaths
    return mc_a - (mc_g - ge_exact)


@pytest.mark.parametrize("ot", ["call", "put"])
def test_matches_control_variate_mc(ot):
    n = 12
    cf = discrete_arithmetic_asian(S, K, T, R, SIG, n_fixings=n, option_type=ot)
    cv = _mc_control_variate(n, ot)
    # Levy is a two-moment approximation; agree to a few tenths of a percent.
    assert cf == pytest.approx(cv, rel=0.01)


def test_cheaper_than_vanilla_call():
    van = price(S, K, T, R, SIG, "call")
    asian = discrete_arithmetic_asian(S, K, T, R, SIG, n_fixings=12, option_type="call")
    assert asian < van


def test_greeks_signs_and_price_field():
    g = discrete_arithmetic_asian_greeks(S, K, T, R, SIG, n_fixings=12,
                                         option_type="call")
    assert g["delta"] > 0.0
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0
    assert g["price"] == pytest.approx(
        discrete_arithmetic_asian(S, K, T, R, SIG, n_fixings=12, option_type="call"),
        abs=1e-9)


def test_validation():
    with pytest.raises(ValueError):
        discrete_arithmetic_asian(S, K, T, R, SIG)
    with pytest.raises(ValueError):
        discrete_arithmetic_asian(S, K, T, R, SIG, fixing_times=[0.5, 2.0])
