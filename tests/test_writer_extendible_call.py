"""Writer-extendible call (Longstaff 1990) — extendible.writer_extendible_call."""

import math

import pytest

from quantforge import writer_extendible_call, writer_extendible_call_greeks
from quantforge.bsm import call_price


S, R, SIG = 100.0, 0.05, 0.25


def _mc(K1, K2, t1, T2, npaths=400000, seed=11):
    import random
    rng = random.Random(seed)
    disc = math.exp(-R * t1)
    acc = 0.0
    for _ in range(npaths):
        z = rng.gauss(0.0, 1.0)
        St1 = S * math.exp((R - 0.5 * SIG * SIG) * t1 + SIG * math.sqrt(t1) * z)
        if St1 > K1:
            acc += St1 - K1
        else:
            acc += call_price(St1, K2, T2 - t1, R, SIG)
    return disc * acc / npaths


@pytest.mark.parametrize("K1,K2,t1,T2", [
    (100.0, 100.0, 0.5, 1.0),
    (100.0, 105.0, 0.5, 1.5),
    (95.0, 100.0, 0.25, 1.0),
])
def test_matches_monte_carlo(K1, K2, t1, T2):
    cf = writer_extendible_call(S, K1, K2, t1, T2, R, SIG)
    mc = _mc(K1, K2, t1, T2)
    assert cf == pytest.approx(mc, rel=0.01)


def test_extension_adds_value_over_plain_call():
    cf = writer_extendible_call(S, 100.0, 100.0, 0.5, 1.0, R, SIG)
    assert cf > call_price(S, 100.0, 0.5, R, SIG)


def test_greeks_price_field_and_delta_fd():
    K1, K2, t1, T2 = 100.0, 100.0, 0.5, 1.0
    g = writer_extendible_call_greeks(S, K1, K2, t1, T2, R, SIG)
    assert g["price"] == pytest.approx(
        writer_extendible_call(S, K1, K2, t1, T2, R, SIG), abs=1e-9)
    h = 1e-4 * S
    fd = (writer_extendible_call(S + h, K1, K2, t1, T2, R, SIG)
          - writer_extendible_call(S - h, K1, K2, t1, T2, R, SIG)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-4)
    assert g["delta"] > 0.0  # call: value rises with spot
    assert g["vega"] > 0.0


def test_validation():
    with pytest.raises(ValueError):
        writer_extendible_call(S, 100.0, 100.0, 1.0, 0.5, R, SIG)  # T2 < t1
