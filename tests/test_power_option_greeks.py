"""Greeks of a power option (power_option_greeks)."""

import pytest

from quantforge import power_option_greeks, power_option, greeks, OptionType


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2


def test_power_one_reduces_to_vanilla_greeks():
    g = power_option_greeks(S, K, T, R, SIG, 1.0, OptionType.CALL)
    van = greeks(S, K, T, R, SIG, OptionType.CALL)
    assert g["delta"] == pytest.approx(van.delta, abs=1e-4)
    assert g["gamma"] == pytest.approx(van.gamma, abs=1e-5)
    assert g["vega"] == pytest.approx(van.vega, abs=1e-2)


def test_delta_matches_finite_difference():
    g = power_option_greeks(S, K, T, R, SIG, 2.0, OptionType.CALL)
    h = 0.01
    fd = (power_option(S + h, K, T, R, SIG, 2.0)
          - power_option(S - h, K, T, R, SIG, 2.0)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, rel=1e-4)


def test_call_greek_signs():
    g = power_option_greeks(S, K, T, R, SIG, 2.0, OptionType.CALL)
    assert g["delta"] > 0.0
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0


def test_put_delta_negative():
    g = power_option_greeks(S, K, T, R, SIG, 2.0, OptionType.PUT)
    assert g["delta"] < 0.0


def test_price_field_matches_power_option():
    g = power_option_greeks(S, K, T, R, SIG, 2.0, OptionType.CALL)
    assert g["price"] == pytest.approx(
        power_option(S, K, T, R, SIG, 2.0, OptionType.CALL), abs=1e-9)


def test_bad_power_raises():
    with pytest.raises(ValueError):
        power_option_greeks(S, K, T, R, SIG, 0.0)
