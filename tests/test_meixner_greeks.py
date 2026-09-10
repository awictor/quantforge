"""Greeks of a Meixner option (meixner_greeks)."""

import pytest

from quantforge import meixner_greeks, meixner_price, OptionType


S, K, T, R = 100.0, 100.0, 1.0, 0.05
A, B, D = 0.3, -0.5, 0.5


def test_delta_matches_finite_difference():
    g = meixner_greeks(S, K, T, R, A, B, D, OptionType.CALL)
    h = 0.01
    fd = (meixner_price(S + h, K, T, R, A, B, D)
          - meixner_price(S - h, K, T, R, A, B, D)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-4)


def test_call_greek_signs():
    g = meixner_greeks(S, K, T, R, A, B, D, OptionType.CALL)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0


def test_put_delta_negative():
    g = meixner_greeks(S, K, T, R, A, B, D, OptionType.PUT)
    assert g["delta"] < 0.0


def test_skew_sensitivity_finite():
    g = meixner_greeks(S, K, T, R, A, B, D, OptionType.CALL)
    assert g["d_b"] == g["d_b"]  # not NaN


def test_price_field_matches_price():
    g = meixner_greeks(S, K, T, R, A, B, D, OptionType.CALL)
    assert g["price"] == pytest.approx(
        meixner_price(S, K, T, R, A, B, D, OptionType.CALL), abs=1e-9)


def test_bad_b_raises():
    with pytest.raises(ValueError):
        meixner_greeks(S, K, T, R, A, 4.0, D)  # b outside (-pi, pi)
