"""Bachelier analytic theta and cash/asset-or-nothing digitals."""

import math

import pytest

from quantforge import (
    bachelier_price, bachelier_theta, bachelier_greeks,
    bachelier_cash_or_nothing, bachelier_asset_or_nothing,
)


F, K, R, SIG, T = 100.0, 100.0, 0.03, 15.0, 1.0


@pytest.mark.parametrize("Fv,Kv,ot", [
    (100.0, 100.0, "call"),
    (105.0, 100.0, "call"),
    (95.0, 100.0, "call"),
    (105.0, 100.0, "put"),
    (95.0, 100.0, "put"),
])
def test_theta_matches_finite_difference(Fv, Kv, ot):
    an = bachelier_theta(Fv, Kv, T, R, SIG, ot)
    ht = 1e-4
    fd = -(bachelier_price(Fv, Kv, T + ht, R, SIG, ot)
           - bachelier_price(Fv, Kv, T - ht, R, SIG, ot)) / (2 * ht)
    assert an == pytest.approx(fd, abs=1e-3)


def test_greeks_uses_analytic_theta():
    g = bachelier_greeks(F, K, T, R, SIG, "call")
    assert g["theta"] == pytest.approx(bachelier_theta(F, K, T, R, SIG, "call"), abs=1e-12)


def test_cash_digital_parity():
    c = bachelier_cash_or_nothing(F, K, T, R, SIG, "call")
    p = bachelier_cash_or_nothing(F, K, T, R, SIG, "put")
    assert c + p == pytest.approx(math.exp(-R * T), abs=1e-12)


def test_cash_digital_matches_monte_carlo():
    import random
    rng = random.Random(3)
    N = 300000
    disc = math.exp(-R * T)
    hits = 0
    for _ in range(N):
        FT = F + SIG * math.sqrt(T) * rng.gauss(0.0, 1.0)
        if FT > K:
            hits += 1
    mc = disc * hits / N
    assert bachelier_cash_or_nothing(F, K, T, R, SIG, "call") == pytest.approx(mc, rel=0.02)


def test_vanilla_decomposition():
    # call = asset-or-nothing - K * cash-or-nothing; put = K*cash - asset.
    c = (bachelier_asset_or_nothing(F, K, T, R, SIG, "call")
         - K * bachelier_cash_or_nothing(F, K, T, R, SIG, "call"))
    p = (K * bachelier_cash_or_nothing(F, K, T, R, SIG, "put")
         - bachelier_asset_or_nothing(F, K, T, R, SIG, "put"))
    assert c == pytest.approx(bachelier_price(F, K, T, R, SIG, "call"), abs=1e-10)
    assert p == pytest.approx(bachelier_price(F, K, T, R, SIG, "put"), abs=1e-10)


def test_asset_digital_matches_monte_carlo():
    import random
    rng = random.Random(5)
    N = 300000
    disc = math.exp(-R * T)
    acc = 0.0
    for _ in range(N):
        FT = F + SIG * math.sqrt(T) * rng.gauss(0.0, 1.0)
        if FT > K:
            acc += FT
    mc = disc * acc / N
    assert bachelier_asset_or_nothing(F, K, T, R, SIG, "call") == pytest.approx(mc, rel=0.02)
