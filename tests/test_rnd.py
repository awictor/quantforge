"""Breeden-Litzenberger risk-neutral density from a smile + payoff pricing."""

import math

import pytest

from quantforge import (
    OptionType,
    risk_neutral_density_from_smile as rnd,
    density_grid_from_smile,
    price_payoff_from_density as pp,
    call_price,
    put_price,
    cash_or_nothing,
)


S0, T, R = 100.0, 1.0, 0.05


def _flat(K):
    return 0.2


def test_density_integrates_to_one():
    ks, d = density_grid_from_smile(S0, T, R, _flat, n=600, width=8.0)
    dK = ks[1] - ks[0]
    mass = sum(max(x, 0.0) for x in d) * dK
    assert mass == pytest.approx(1.0, abs=5e-3)


def test_reprices_call():
    c = pp(S0, T, R, _flat, lambda K: max(K - 100, 0.0), n=800, width=8.0)
    assert c == pytest.approx(call_price(S0, 100, T, R, 0.2), abs=1e-2)


def test_reprices_put():
    p = pp(S0, T, R, _flat, lambda K: max(100 - K, 0.0), n=800, width=8.0)
    assert p == pytest.approx(put_price(S0, 100, T, R, 0.2), abs=1e-2)


def test_reprices_digital():
    dig = pp(S0, T, R, _flat, lambda K: 1.0 if K > 100 else 0.0,
             n=800, width=8.0)
    assert dig == pytest.approx(
        cash_or_nothing(S0, 100, T, R, 0.2, OptionType.CALL), abs=1e-2)


def test_forward_recovers_spot():
    # payoff(K) = K prices the discounted expected terminal spot = S0 (q=0).
    fwd = pp(S0, T, R, _flat, lambda K: K, n=800, width=8.0)
    assert fwd == pytest.approx(S0, abs=1e-1)


def test_atm_density_positive():
    assert rnd(S0, T, R, _flat, S0) > 0.0


def test_arbitrage_smile_gives_negative_density():
    # A concave dip in the smile creates butterfly arbitrage -> negative density.
    def arb(K):
        return max(0.01, 0.2 - 2.0 * abs(math.log(K / 100)))
    # Somewhere in the wing the density goes negative.
    dens = [rnd(S0, T, R, arb, K) for K in (85, 90, 110, 115)]
    assert min(dens) < 0.0
