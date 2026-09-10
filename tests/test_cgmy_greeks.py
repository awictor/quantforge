"""Greeks of a CGMY option (cgmy_greeks)."""

import pytest

from quantforge import cgmy_greeks, cgmy_price, OptionType


S, K, T, R = 100.0, 100.0, 1.0, 0.05
C, G, M, Y = 0.1, 5.0, 5.0, 0.5


def test_delta_matches_finite_difference():
    g = cgmy_greeks(S, K, T, R, C, G, M, Y, OptionType.CALL)
    h = 0.01
    fd = (cgmy_price(S + h, K, T, R, C, G, M, Y)
          - cgmy_price(S - h, K, T, R, C, G, M, Y)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-4)


def test_call_greek_signs():
    g = cgmy_greeks(S, K, T, R, C, G, M, Y, OptionType.CALL)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0


def test_put_delta_negative():
    g = cgmy_greeks(S, K, T, R, C, G, M, Y, OptionType.PUT)
    assert g["delta"] < 0.0


def test_tail_sensitivity_finite():
    g = cgmy_greeks(S, K, T, R, C, G, M, Y, OptionType.CALL)
    assert g["d_Y"] == g["d_Y"]  # not NaN


def test_price_field_matches_price():
    g = cgmy_greeks(S, K, T, R, C, G, M, Y, OptionType.CALL)
    assert g["price"] == pytest.approx(
        cgmy_price(S, K, T, R, C, G, M, Y, OptionType.CALL), abs=1e-9)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        cgmy_greeks(S, K, T, R, C, G, M, 2.5)  # Y >= 2
