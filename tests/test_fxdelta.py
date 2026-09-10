"""FX delta-space quoting conventions."""

import math

import pytest

from quantforge import (
    atm_dns_strike,
    strike_from_delta,
    delta_from_strike,
    rr_bf_to_pillars,
)


F, T, SIGMA = 1.30, 1.0, 0.10


@pytest.mark.parametrize("delta,is_call", [
    (0.25, True), (-0.25, False), (0.10, True), (-0.10, False),
])
def test_forward_delta_round_trip(delta, is_call):
    K = strike_from_delta(F, T, SIGMA, delta, is_call)
    got = delta_from_strike(F, T, SIGMA, K, is_call)
    assert got == pytest.approx(delta, abs=1e-10)


def test_premium_adjusted_round_trip():
    K = strike_from_delta(F, T, SIGMA, 0.25, True, premium_adjusted=True)
    got = delta_from_strike(F, T, SIGMA, K, True, premium_adjusted=True)
    assert got == pytest.approx(0.25, abs=1e-8)


def test_spot_delta_round_trip():
    K = strike_from_delta(F, T, SIGMA, 0.25, True, spot_delta=True, r_for=0.02)
    got = delta_from_strike(F, T, SIGMA, K, True, spot_delta=True, r_for=0.02)
    assert got == pytest.approx(0.25, abs=1e-8)


def test_atm_dns_is_delta_neutral_straddle():
    K = atm_dns_strike(F, T, SIGMA)
    dc = delta_from_strike(F, T, SIGMA, K, True)
    dp = delta_from_strike(F, T, SIGMA, K, False)
    assert (dc + dp) == pytest.approx(0.0, abs=1e-10)


def test_rr_bf_pillars_reproduce_quotes():
    atm, rr, bf = 0.10, -0.015, 0.004
    Kp, sp, Ka, sa, Kc, sc = rr_bf_to_pillars(F, T, atm, rr, bf, 0.25)
    assert (sc - sp) == pytest.approx(rr, abs=1e-12)          # risk reversal
    assert ((sc + sp) / 2 - sa) == pytest.approx(bf, abs=1e-12)  # butterfly
    assert Kp < Ka < Kc                                       # strike ordering


def test_negative_rr_puts_downside_vol_higher():
    # A negative 25d risk reversal makes the 25d put vol exceed the call vol.
    _, sp, _, _, _, sc = rr_bf_to_pillars(F, T, 0.10, -0.02, 0.005, 0.25)
    assert sp > sc


def test_atm_strike_above_forward_for_positive_vol():
    assert atm_dns_strike(F, T, SIGMA) > F
    assert atm_dns_strike(F, T, 0.0) == pytest.approx(F)
