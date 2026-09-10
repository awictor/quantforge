"""Greeks of a Kou double-exponential jump-diffusion option (kou_greeks)."""

import pytest

from quantforge import kou_greeks, kou_price, greeks, OptionType


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2
LAM, P, ETA1, ETA2 = 1.0, 0.4, 10.0, 5.0


def test_no_jumps_reduces_to_vanilla():
    g = kou_greeks(S, K, T, R, SIG, 0.0, 0.5, 10.0, 10.0, OptionType.CALL)
    van = greeks(S, K, T, R, SIG, OptionType.CALL)
    assert g["delta"] == pytest.approx(van.delta, abs=1e-4)
    assert g["gamma"] == pytest.approx(van.gamma, abs=1e-5)
    assert g["vega"] == pytest.approx(van.vega, abs=1e-2)


def test_delta_matches_finite_difference():
    g = kou_greeks(S, K, T, R, SIG, LAM, P, ETA1, ETA2, OptionType.CALL)
    h = 0.01
    fd = (kou_price(S + h, K, T, R, SIG, LAM, P, ETA1, ETA2)
          - kou_price(S - h, K, T, R, SIG, LAM, P, ETA1, ETA2)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-4)


def test_call_greek_signs():
    g = kou_greeks(S, K, T, R, SIG, LAM, P, ETA1, ETA2, OptionType.CALL)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0


def test_put_delta_negative():
    g = kou_greeks(S, K, T, R, SIG, LAM, P, ETA1, ETA2, OptionType.PUT)
    assert g["delta"] < 0.0


def test_price_field_matches_price():
    g = kou_greeks(S, K, T, R, SIG, LAM, P, ETA1, ETA2, OptionType.CALL)
    assert g["price"] == pytest.approx(
        kou_price(S, K, T, R, SIG, LAM, P, ETA1, ETA2, OptionType.CALL),
        abs=1e-9)


def test_bad_lambda_raises():
    with pytest.raises(ValueError):
        kou_greeks(S, K, T, R, SIG, -1.0, P, ETA1, ETA2)
