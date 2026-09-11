"""BSM implied-vol smile produced by the CEV model (cev.cev_smile)."""

import math

import pytest

from quantforge import cev_smile, cev_price, implied_volatility


S, T, R, SIG = 100.0, 1.0, 0.05, 0.2
STRIKES = [80.0, 90.0, 100.0, 110.0, 120.0]


def test_atm_vol_near_sigma():
    # sigma is calibrated to the ATM instantaneous vol, so the ATM implied vol
    # sits close to it.
    sm = cev_smile(S, [S * math.exp(R * T)], T, R, SIG, 0.5)
    assert sm[0][1] == pytest.approx(SIG, abs=5e-3)


def test_downward_skew():
    sm = dict(cev_smile(S, [80.0, 120.0], T, R, SIG, 0.3))
    # log-moneyness keys: left wing (K=80) < 0, right wing (K=120) > 0.
    left = [v for k, v in sm.items() if k < 0][0]
    right = [v for k, v in sm.items() if k > 0][0]
    assert left > right  # CEV beta < 1 gives a downward skew


def test_lower_beta_steeper_skew():
    def slope(beta):
        sm = dict(cev_smile(S, [80.0, 120.0], T, R, SIG, beta))
        ks = sorted(sm)
        return (sm[ks[-1]] - sm[ks[0]]) / (ks[-1] - ks[0])
    # More negative slope for smaller beta.
    assert slope(0.2) < slope(0.5) < slope(0.9)


def test_output_sorted_and_paired():
    sm = cev_smile(S, list(reversed(STRIKES)), T, R, SIG, 0.5)
    ks = [k for k, _ in sm]
    assert ks == sorted(ks)
    assert len(sm) == len(STRIKES)


def test_reprices_to_cev():
    # Each smile vol, put back through Black-Scholes, must reproduce the CEV
    # price it came from.
    beta = 0.4
    F = S * math.exp(R * T)
    for lm, iv in cev_smile(S, STRIKES, T, R, SIG, beta):
        K = F * math.exp(lm)
        cev = cev_price(S, K, T, R, SIG, beta, "call")
        recov = implied_volatility(cev, S, K, T, R, "call", b=R)
        assert iv == pytest.approx(recov, abs=1e-9)
