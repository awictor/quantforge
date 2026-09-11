"""Average-strike (floating-strike) discrete geometric Asian (exchange option)."""

import math

import pytest

from quantforge import (
    average_strike_geometric_asian, average_strike_geometric_asian_greeks,
)
from quantforge.bsm import OptionType


S, R, SIG = 100.0, 0.05, 0.2
T, N = 1.0, 12


def _E_ST_E_G(S, t, r, sig, n):
    times = [t * i / n for i in range(1, n + 1)]
    mean_t = sum(times) / n
    m_G = math.log(S) + (r - 0.5 * sig * sig) * mean_t
    dbl = sum((ti if ti < tj else tj) for ti in times for tj in times)
    v_G = sig * sig * dbl / (n * n)
    return S * math.exp(r * t), math.exp(m_G + 0.5 * v_G)


def test_put_call_parity_is_exchange_forward():
    c = average_strike_geometric_asian(S, T, R, SIG, N, option_type="call")
    p = average_strike_geometric_asian(S, T, R, SIG, N, option_type="put")
    E_ST, E_G = _E_ST_E_G(S, T, R, SIG, N)
    assert c - p == pytest.approx(math.exp(-R * T) * (E_ST - E_G), abs=1e-9)


def _mc(ot, npaths=1000000, seed=21):
    import random
    rng = random.Random(seed)
    disc = math.exp(-R * T)
    acc = 0.0
    for _ in range(npaths):
        z = [rng.gauss(0.0, 1.0) for _ in range(N)]
        for sgn in (1.0, -1.0):
            logS = math.log(S)
            prevt = 0.0
            gsum = 0.0
            for i in range(1, N + 1):
                ti = T * i / N
                dt = ti - prevt
                prevt = ti
                logS += (R - 0.5 * SIG * SIG) * dt + SIG * math.sqrt(dt) * sgn * z[i - 1]
                gsum += logS
            ST = math.exp(logS)
            G = math.exp(gsum / N)
            pay = max(ST - G, 0.0) if ot == "call" else max(G - ST, 0.0)
            acc += pay
    return disc * acc / (2 * npaths)


@pytest.mark.slow
@pytest.mark.parametrize("ot", ["call", "put"])
def test_matches_monte_carlo(ot):
    cf = average_strike_geometric_asian(S, T, R, SIG, N, option_type=ot)
    mc = _mc(ot)
    assert cf == pytest.approx(mc, rel=0.01)


def test_explicit_times_equal_equally_spaced():
    times = [T * i / N for i in range(1, N + 1)]
    a = average_strike_geometric_asian(S, T, R, SIG, N, option_type="call")
    b = average_strike_geometric_asian(S, T, R, SIG, fixing_times=times,
                                       option_type="call")
    assert a == pytest.approx(b, abs=1e-12)


def test_positive_prices():
    assert average_strike_geometric_asian(S, T, R, SIG, N, option_type="call") > 0.0
    assert average_strike_geometric_asian(S, T, R, SIG, N, option_type="put") > 0.0


def test_greeks_price_field_and_gamma_sign():
    g = average_strike_geometric_asian_greeks(S, T, R, SIG, N, option_type="call")
    assert g["price"] == pytest.approx(
        average_strike_geometric_asian(S, T, R, SIG, N, option_type="call"), abs=1e-9)
    assert g["gamma"] > 0.0   # convex payoff
    assert g["vega"] > 0.0    # more dispersion -> more value


def test_validation():
    with pytest.raises(ValueError):
        average_strike_geometric_asian(S, T, R, SIG)  # neither n nor times
    with pytest.raises(ValueError):
        average_strike_geometric_asian(S, T, R, SIG, fixing_times=[0.5, 1.5])  # > t
