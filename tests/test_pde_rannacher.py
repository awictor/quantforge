"""Rannacher time-stepping in the Crank-Nicolson PDE solver."""

import pytest

from quantforge import (
    OptionType,
    crank_nicolson_price,
    crank_nicolson_greeks,
    call_price,
    gamma as bs_gamma,
)


S, K, T, R, SIGMA = 100.0, 100.0, 1.0, 0.05, 0.2


def test_rannacher_price_matches_black_scholes():
    p = crank_nicolson_price(S, K, T, R, SIGMA, OptionType.CALL,
                             n_space=300, n_time=300, rannacher=2)
    assert p == pytest.approx(call_price(S, K, T, R, SIGMA), abs=1e-2)


def test_rannacher_zero_recovers_pure_cn():
    # rannacher=0 is pure Crank-Nicolson; both still match BS.
    p_cn = crank_nicolson_price(S, K, T, R, SIGMA, n_space=300, n_time=300,
                                rannacher=0)
    p_r = crank_nicolson_price(S, K, T, R, SIGMA, n_space=300, n_time=300,
                               rannacher=4)
    bs = call_price(S, K, T, R, SIGMA)
    assert p_cn == pytest.approx(bs, abs=1e-2)
    assert p_r == pytest.approx(bs, abs=1e-2)


def test_rannacher_gamma_no_worse_at_coarse_time_grid():
    # With few time steps the payoff kink makes pure-CN gamma oscillate;
    # Rannacher damping should not be worse than pure CN near the strike.
    ref = bs_gamma(S, K, 0.25, R, SIGMA)
    g_cn = crank_nicolson_greeks(S, K, 0.25, R, SIGMA, n_space=200, n_time=25,
                                 rannacher=0)["gamma"]
    g_r = crank_nicolson_greeks(S, K, 0.25, R, SIGMA, n_space=200, n_time=25,
                                rannacher=2)["gamma"]
    assert abs(g_r - ref) <= abs(g_cn - ref) + 1e-6


def test_rannacher_greeks_still_match_black_scholes():
    g = crank_nicolson_greeks(S, K, T, R, SIGMA, OptionType.CALL,
                              n_space=400, n_time=400, rannacher=2)
    assert g["gamma"] == pytest.approx(bs_gamma(S, K, T, R, SIGMA), abs=1e-3)
