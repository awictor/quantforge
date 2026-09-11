"""Bjerksund-Stensland (2014) spread-option approximation (spread_option_bs)."""

import math

import pytest

from quantforge import (
    spread_option_bs, spread_option, exchange_option, OptionType,
)
from quantforge.montecarlo import spread_option_lhs_mc


T, R = 1.0, 0.03


def test_zero_strike_equals_margrabe():
    # At K = 0 the payoff is max(S1 - S2, 0): the exact Margrabe value.
    S1, S2, s1, s2, rho = 100.0, 100.0, 0.25, 0.25, 0.7
    bs = spread_option_bs(S1, S2, 0.0, T, R, s1, s2, rho)
    marg = exchange_option(S1, S2, T, s1, s2, rho)
    assert bs == pytest.approx(marg, abs=1e-9)


def test_put_call_parity():
    # C - P = disc*(F1 - F2 - K).
    S1, S2, K, s1, s2, rho = 100.0, 90.0, 10.0, 0.3, 0.35, 0.3
    F1 = S1 * math.exp(R * T)
    F2 = S2 * math.exp(R * T)
    disc = math.exp(-R * T)
    c = spread_option_bs(S1, S2, K, T, R, s1, s2, rho, option_type=OptionType.CALL)
    p = spread_option_bs(S1, S2, K, T, R, s1, s2, rho, option_type=OptionType.PUT)
    assert c - p == pytest.approx(disc * (F1 - F2 - K), abs=1e-9)


@pytest.mark.parametrize("S1,S2,K,s1,s2,rho", [
    (100.0, 96.0, 4.0, 0.2, 0.3, 0.5),
    (100.0, 90.0, 10.0, 0.3, 0.35, 0.3),
    (120.0, 100.0, 15.0, 0.4, 0.2, -0.2),
])
def test_matches_monte_carlo(S1, S2, K, s1, s2, rho):
    bs = spread_option_bs(S1, S2, K, T, R, s1, s2, rho)
    mc = spread_option_lhs_mc(S1, S2, K, T, R, s1, s2, rho,
                              n_paths=200000, seed=7)
    assert bs == pytest.approx(mc.price, abs=4.0 * mc.std_error + 1e-3)


def test_accuracy_versus_kirk_at_wide_strike():
    # Bjerksund-Stensland should sit at least as close to MC as Kirk at a wide,
    # dispersed-vol strike where Kirk's blended-vol assumption strains.
    S1, S2, K, s1, s2, rho = 120.0, 100.0, 15.0, 0.4, 0.2, -0.2
    mc = spread_option_lhs_mc(S1, S2, K, T, R, s1, s2, rho,
                              n_paths=400000, seed=11)
    bs = spread_option_bs(S1, S2, K, T, R, s1, s2, rho)
    kirk = spread_option(S1, S2, K, T, R, s1, s2, rho)
    assert abs(bs - mc.price) <= abs(kirk - mc.price) + 2.0 * mc.std_error


def test_call_decreasing_in_correlation():
    # Higher correlation shrinks the spread vol, lowering the call.
    S1, S2, K, s1, s2 = 100.0, 95.0, 5.0, 0.3, 0.3
    lo = spread_option_bs(S1, S2, K, T, R, s1, s2, -0.5)
    hi = spread_option_bs(S1, S2, K, T, R, s1, s2, 0.9)
    assert hi < lo


def test_positive_prices_negative_strike():
    # Negative K (S2 gets a subsidy) is allowed as long as F2 + K > 0.
    v = spread_option_bs(100.0, 100.0, -20.0, T, R, 0.25, 0.25, 0.4)
    assert v > 0.0
