"""Greeks of a Merton jump-diffusion option (merton_jump_greeks)."""

import pytest

from quantforge import merton_jump_greeks, merton_jump_price, greeks, OptionType


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2
LAM, MU_J, SIG_J = 0.5, -0.1, 0.15


def test_no_jumps_reduces_to_vanilla():
    g = merton_jump_greeks(S, K, T, R, SIG, 0.0, 0.0, 0.0, OptionType.CALL)
    van = greeks(S, K, T, R, SIG, OptionType.CALL)
    assert g["delta"] == pytest.approx(van.delta, abs=1e-4)
    assert g["gamma"] == pytest.approx(van.gamma, abs=1e-5)
    assert g["vega"] == pytest.approx(van.vega, abs=1e-2)


def test_delta_matches_finite_difference():
    g = merton_jump_greeks(S, K, T, R, SIG, LAM, MU_J, SIG_J, OptionType.CALL)
    h = 0.01
    fd = (merton_jump_price(S + h, K, T, R, SIG, LAM, MU_J, SIG_J)
          - merton_jump_price(S - h, K, T, R, SIG, LAM, MU_J, SIG_J)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-5)


def test_call_greek_signs():
    g = merton_jump_greeks(S, K, T, R, SIG, LAM, MU_J, SIG_J, OptionType.CALL)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0


def test_put_delta_negative():
    g = merton_jump_greeks(S, K, T, R, SIG, LAM, MU_J, SIG_J, OptionType.PUT)
    assert g["delta"] < 0.0


def test_price_field_matches_price():
    g = merton_jump_greeks(S, K, T, R, SIG, LAM, MU_J, SIG_J, OptionType.CALL)
    assert g["price"] == pytest.approx(
        merton_jump_price(S, K, T, R, SIG, LAM, MU_J, SIG_J, OptionType.CALL),
        abs=1e-12)


def test_bad_lambda_raises():
    with pytest.raises(ValueError):
        merton_jump_greeks(S, K, T, R, SIG, -1.0, MU_J, SIG_J)
