"""Crank-Nicolson PDE Greeks read off the grid vs Black-Scholes."""

import pytest

from quantforge import (
    OptionType,
    crank_nicolson_greeks,
    crank_nicolson_price,
    delta as bs_delta,
    gamma as bs_gamma,
    theta as bs_theta,
)


S, K, T, R, SIGMA = 100.0, 100.0, 1.0, 0.05, 0.2


def test_call_greeks_match_black_scholes():
    g = crank_nicolson_greeks(S, K, T, R, SIGMA, OptionType.CALL,
                              n_space=400, n_time=400)
    assert g["delta"] == pytest.approx(bs_delta(S, K, T, R, SIGMA), abs=2e-3)
    assert g["gamma"] == pytest.approx(bs_gamma(S, K, T, R, SIGMA), abs=1e-3)
    assert g["theta"] == pytest.approx(bs_theta(S, K, T, R, SIGMA), abs=5e-2)


def test_put_delta_matches_black_scholes():
    g = crank_nicolson_greeks(S, K, T, R, SIGMA, OptionType.PUT,
                              n_space=400, n_time=400)
    assert g["delta"] == pytest.approx(
        bs_delta(S, K, T, R, SIGMA, OptionType.PUT), abs=2e-3)
    assert g["gamma"] == pytest.approx(bs_gamma(S, K, T, R, SIGMA), abs=1e-3)


def test_price_field_matches_price_function():
    g = crank_nicolson_greeks(S, K, T, R, SIGMA, OptionType.CALL,
                              n_space=300, n_time=300)
    p = crank_nicolson_price(S, K, T, R, SIGMA, OptionType.CALL,
                             n_space=300, n_time=300)
    assert g["price"] == pytest.approx(p, abs=1e-9)


def test_call_delta_in_unit_interval():
    g = crank_nicolson_greeks(S, K, T, R, SIGMA, OptionType.CALL,
                              n_space=200, n_time=200)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0
    assert g["theta"] < 0.0


@pytest.mark.slow
def test_american_put_greeks_reasonable():
    g = crank_nicolson_greeks(S, K, T, R, SIGMA, OptionType.PUT, american=True,
                              n_space=300, n_time=300)
    # American put: negative delta, positive gamma.
    assert -1.0 < g["delta"] < 0.0
    assert g["gamma"] > 0.0


def test_requires_positive_time():
    with pytest.raises(ValueError):
        crank_nicolson_greeks(S, K, 0.0, R, SIGMA)
