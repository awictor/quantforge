"""Holder-extendible put (Longstaff 1990) — extendible.holder_extendible_put."""

import math

import pytest

from quantforge import holder_extendible_put, holder_extendible_put_greeks
from quantforge.bsm import put_price


S, R, SIG = 100.0, 0.05, 0.25


def _mc(K1, K2, t1, T2, A, npaths=400000, seed=13):
    import random
    rng = random.Random(seed)
    disc = math.exp(-R * t1)
    acc = 0.0
    for _ in range(npaths):
        z = rng.gauss(0.0, 1.0)
        St1 = S * math.exp((R - 0.5 * SIG * SIG) * t1 + SIG * math.sqrt(t1) * z)
        acc += max(K1 - St1, put_price(St1, K2, T2 - t1, R, SIG) - A, 0.0)
    return disc * acc / npaths


@pytest.mark.parametrize("K1,K2,t1,T2,A", [
    (100.0, 100.0, 0.5, 1.0, 1.0),
    (100.0, 95.0, 0.5, 1.5, 2.0),
    (105.0, 100.0, 0.25, 1.0, 0.5),
])
def test_matches_monte_carlo(K1, K2, t1, T2, A):
    cf = holder_extendible_put(S, K1, K2, t1, T2, R, SIG, A)
    mc = _mc(K1, K2, t1, T2, A)
    assert cf == pytest.approx(mc, rel=0.01)


def test_huge_fee_reduces_to_vanilla_put():
    cf = holder_extendible_put(S, 100.0, 100.0, 0.5, 1.0, R, SIG, 1e6)
    assert cf == pytest.approx(put_price(S, 100.0, 0.5, R, SIG), abs=1e-6)


def test_extension_adds_value():
    cf = holder_extendible_put(S, 100.0, 100.0, 0.5, 1.0, R, SIG, 0.5)
    assert cf > put_price(S, 100.0, 0.5, R, SIG)


def test_lower_fee_worth_more():
    cheap = holder_extendible_put(S, 100.0, 100.0, 0.5, 1.0, R, SIG, 0.5)
    dear = holder_extendible_put(S, 100.0, 100.0, 0.5, 1.0, R, SIG, 3.0)
    assert cheap > dear


def test_greeks_price_field_and_delta_fd():
    K1, K2, t1, T2, A = 100.0, 100.0, 0.5, 1.0, 1.0
    g = holder_extendible_put_greeks(S, K1, K2, t1, T2, R, SIG, A)
    assert g["price"] == pytest.approx(
        holder_extendible_put(S, K1, K2, t1, T2, R, SIG, A), abs=1e-9)
    h = 1e-4 * S
    fd = (holder_extendible_put(S + h, K1, K2, t1, T2, R, SIG, A)
          - holder_extendible_put(S - h, K1, K2, t1, T2, R, SIG, A)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-4)
    assert g["delta"] < 0.0  # put falls as spot rises
    assert g["vega"] > 0.0


def test_validation():
    with pytest.raises(ValueError):
        holder_extendible_put(S, 100.0, 100.0, 1.0, 0.5, R, SIG, 1.0)  # T2 < t1
    with pytest.raises(ValueError):
        holder_extendible_put(S, 100.0, 100.0, 0.5, 1.0, R, SIG, -1.0)  # A < 0
