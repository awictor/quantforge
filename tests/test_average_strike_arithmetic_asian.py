"""Average-strike (floating-strike) discrete arithmetic Asian (exchange option)."""

import math

import pytest

from quantforge import (
    average_strike_arithmetic_asian, average_strike_arithmetic_asian_greeks,
    average_strike_geometric_asian,
)


S, R, SIG = 100.0, 0.05, 0.2
T, N = 1.0, 12


def test_put_call_parity_is_exchange_forward():
    c = average_strike_arithmetic_asian(S, T, R, SIG, N, option_type="call")
    p = average_strike_arithmetic_asian(S, T, R, SIG, N, option_type="put")
    E_ST = S * math.exp(R * T)
    E_A = (S / N) * sum(math.exp(R * T * i / N) for i in range(1, N + 1))
    assert c - p == pytest.approx(math.exp(-R * T) * (E_ST - E_A), abs=1e-9)


def test_arithmetic_call_cheaper_than_geometric():
    # A >= G (AM-GM) => arithmetic strike is higher => average-strike call
    # is worth less than the geometric-strike call.
    ar = average_strike_arithmetic_asian(S, T, R, SIG, N, option_type="call")
    ge = average_strike_geometric_asian(S, T, R, SIG, N, option_type="call")
    assert ar < ge


def _mc(ot, npaths=1000000, seed=31):
    import random
    rng = random.Random(seed)
    disc = math.exp(-R * T)
    acc = 0.0
    for _ in range(npaths):
        z = [rng.gauss(0.0, 1.0) for _ in range(N)]
        for sgn in (1.0, -1.0):
            logS = math.log(S)
            prevt = 0.0
            asum = 0.0
            for i in range(1, N + 1):
                ti = T * i / N
                dt = ti - prevt
                prevt = ti
                logS += (R - 0.5 * SIG * SIG) * dt + SIG * math.sqrt(dt) * sgn * z[i - 1]
                asum += math.exp(logS)
            ST = math.exp(logS)
            A = asum / N
            pay = max(ST - A, 0.0) if ot == "call" else max(A - ST, 0.0)
            acc += pay
    return disc * acc / (2 * npaths)


@pytest.mark.slow
@pytest.mark.parametrize("ot", ["call", "put"])
def test_matches_monte_carlo(ot):
    # Levy two-moment match: agree to ~1% (a small stable approximation bias).
    cf = average_strike_arithmetic_asian(S, T, R, SIG, N, option_type=ot)
    mc = _mc(ot)
    assert cf == pytest.approx(mc, rel=0.02)


def test_explicit_times_equal_equally_spaced():
    times = [T * i / N for i in range(1, N + 1)]
    a = average_strike_arithmetic_asian(S, T, R, SIG, N, option_type="call")
    b = average_strike_arithmetic_asian(S, T, R, SIG, fixing_times=times,
                                        option_type="call")
    assert a == pytest.approx(b, abs=1e-12)


def test_greeks_price_field_and_homogeneity():
    # The payoff is homogeneous of degree 1 in S (both S_T and A scale with S),
    # so the price is linear in S: gamma is exactly zero and delta = price / S.
    g = average_strike_arithmetic_asian_greeks(S, T, R, SIG, N, option_type="call")
    price = average_strike_arithmetic_asian(S, T, R, SIG, N, option_type="call")
    assert g["price"] == pytest.approx(price, abs=1e-9)
    assert g["gamma"] == pytest.approx(0.0, abs=1e-6)
    assert g["delta"] == pytest.approx(price / S, rel=1e-4)
    assert g["vega"] > 0.0


def test_validation():
    with pytest.raises(ValueError):
        average_strike_arithmetic_asian(S, T, R, SIG)
    with pytest.raises(ValueError):
        average_strike_arithmetic_asian(S, T, R, SIG, fixing_times=[0.5, 2.0])
