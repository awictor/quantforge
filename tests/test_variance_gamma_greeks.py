"""Greeks of a Variance-Gamma option (variance_gamma_greeks)."""

import pytest

from quantforge import (
    variance_gamma_greeks,
    variance_gamma_price,
    greeks,
    OptionType,
)


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2
NU, THETA = 0.2, -0.1


def test_small_nu_near_black_scholes():
    # As nu -> 0 the VG price approaches Black-Scholes; a small nu gives a delta
    # close to the BSM delta.
    g = variance_gamma_greeks(S, K, T, R, SIG, 0.02, -0.1, OptionType.CALL)
    van = greeks(S, K, T, R, SIG, OptionType.CALL)
    assert g["delta"] == pytest.approx(van.delta, abs=0.01)


def test_delta_matches_finite_difference():
    g = variance_gamma_greeks(S, K, T, R, SIG, NU, THETA, OptionType.CALL)
    h = 0.01
    fd = (variance_gamma_price(S + h, K, T, R, SIG, NU, THETA)
          - variance_gamma_price(S - h, K, T, R, SIG, NU, THETA)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-4)


def test_call_greek_signs():
    g = variance_gamma_greeks(S, K, T, R, SIG, NU, THETA, OptionType.CALL)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0


def test_put_delta_negative():
    g = variance_gamma_greeks(S, K, T, R, SIG, NU, THETA, OptionType.PUT)
    assert g["delta"] < 0.0


def test_price_field_matches_price():
    g = variance_gamma_greeks(S, K, T, R, SIG, NU, THETA, OptionType.CALL)
    assert g["price"] == pytest.approx(
        variance_gamma_price(S, K, T, R, SIG, NU, THETA, OptionType.CALL),
        abs=1e-9)


def test_bad_nu_raises():
    with pytest.raises(ValueError):
        variance_gamma_greeks(S, K, T, R, SIG, 0.0, THETA)
