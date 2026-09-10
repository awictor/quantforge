"""Bundled Bachelier Greeks (bachelier_greeks)."""

import pytest

from quantforge import (
    bachelier_greeks,
    bachelier_price,
    bachelier_delta,
    bachelier_gamma,
    bachelier_vega,
    OptionType,
)


F, K, T, R, SIG = 100.0, 100.0, 1.0, 0.02, 15.0


def test_delta_gamma_vega_match_closed_forms():
    g = bachelier_greeks(F, K, T, R, SIG, OptionType.CALL)
    assert g["delta"] == pytest.approx(
        bachelier_delta(F, K, T, R, SIG, OptionType.CALL), abs=1e-12)
    assert g["gamma"] == pytest.approx(bachelier_gamma(F, K, T, R, SIG), abs=1e-12)
    assert g["vega"] == pytest.approx(bachelier_vega(F, K, T, R, SIG), abs=1e-12)


def test_theta_matches_finite_difference():
    g = bachelier_greeks(F, K, T, R, SIG, OptionType.CALL)
    h = 1e-4
    fd = -(bachelier_price(F, K, T + h, R, SIG)
           - bachelier_price(F, K, T - h, R, SIG)) / (2 * h)
    assert g["theta"] == pytest.approx(fd, abs=1e-4)


def test_atm_call_delta_is_half_at_zero_rate():
    g = bachelier_greeks(F, K, T, 0.0, SIG, OptionType.CALL)
    assert g["delta"] == pytest.approx(0.5, abs=1e-9)


def test_put_delta_negative():
    g = bachelier_greeks(F, K, T, R, SIG, OptionType.PUT)
    assert g["delta"] < 0.0


def test_price_field_matches_price():
    g = bachelier_greeks(F, K, T, R, SIG, OptionType.CALL)
    assert g["price"] == pytest.approx(
        bachelier_price(F, K, T, R, SIG, OptionType.CALL), abs=1e-12)
