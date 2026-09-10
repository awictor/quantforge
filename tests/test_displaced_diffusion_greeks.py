"""Greeks of a displaced-diffusion option (displaced_diffusion_greeks)."""

import pytest

from quantforge import (
    displaced_diffusion_greeks,
    displaced_diffusion_price,
    greeks,
    OptionType,
)


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2


def test_zero_shift_reduces_to_vanilla():
    g = displaced_diffusion_greeks(S, K, T, R, SIG, 0.0, OptionType.CALL)
    van = greeks(S, K, T, R, SIG, OptionType.CALL)
    assert g["delta"] == pytest.approx(van.delta, abs=1e-4)
    assert g["gamma"] == pytest.approx(van.gamma, abs=1e-5)
    assert g["vega"] == pytest.approx(van.vega, abs=1e-2)


def test_delta_matches_finite_difference():
    g = displaced_diffusion_greeks(S, K, T, R, SIG, 50.0, OptionType.CALL)
    h = 0.01
    fd = (displaced_diffusion_price(S + h, K, T, R, SIG, 50.0)
          - displaced_diffusion_price(S - h, K, T, R, SIG, 50.0)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-5)


def test_call_greek_signs():
    g = displaced_diffusion_greeks(S, K, T, R, SIG, 50.0, OptionType.CALL)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0


def test_put_delta_negative():
    g = displaced_diffusion_greeks(S, K, T, R, SIG, 50.0, OptionType.PUT)
    assert g["delta"] < 0.0


def test_price_field_matches_price():
    g = displaced_diffusion_greeks(S, K, T, R, SIG, 50.0, OptionType.CALL)
    assert g["price"] == pytest.approx(
        displaced_diffusion_price(S, K, T, R, SIG, 50.0, OptionType.CALL),
        abs=1e-12)


def test_bad_shift_raises():
    with pytest.raises(ValueError):
        displaced_diffusion_greeks(S, K, T, R, SIG, -200.0)
