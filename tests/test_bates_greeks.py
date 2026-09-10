"""Greeks of a Bates (Heston + jumps) option (bates_greeks)."""

import pytest

from quantforge import bates_greeks, bates_price, heston_price, OptionType


S, K, T, R = 100.0, 100.0, 1.0, 0.03
V0, KAPPA, THETA, XI, RHO = 0.04, 1.5, 0.04, 0.5, -0.7
LAM, MU_J, SIG_J = 0.5, -0.1, 0.15


def test_no_jumps_delta_matches_heston():
    g = bates_greeks(S, K, T, R, V0, KAPPA, THETA, XI, RHO, 0.0, 0.0, 0.0,
                     OptionType.CALL)
    h = 0.5
    hd = (heston_price(S + h, K, T, R, V0, KAPPA, THETA, XI, RHO, OptionType.CALL)
          - heston_price(S - h, K, T, R, V0, KAPPA, THETA, XI, RHO,
                         OptionType.CALL)) / (2 * h)
    assert g["delta"] == pytest.approx(hd, abs=1e-3)


def test_delta_matches_finite_difference():
    g = bates_greeks(S, K, T, R, V0, KAPPA, THETA, XI, RHO, LAM, MU_J, SIG_J,
                     OptionType.CALL)
    h = 0.5
    fd = (bates_price(S + h, K, T, R, V0, KAPPA, THETA, XI, RHO, LAM, MU_J, SIG_J)
          - bates_price(S - h, K, T, R, V0, KAPPA, THETA, XI, RHO, LAM, MU_J,
                        SIG_J)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-3)


def test_call_signs():
    g = bates_greeks(S, K, T, R, V0, KAPPA, THETA, XI, RHO, LAM, MU_J, SIG_J,
                     OptionType.CALL)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0
    assert g["vega_v0"] > 0.0


def test_price_field_matches_price():
    g = bates_greeks(S, K, T, R, V0, KAPPA, THETA, XI, RHO, LAM, MU_J, SIG_J,
                     OptionType.CALL)
    assert g["price"] == pytest.approx(
        bates_price(S, K, T, R, V0, KAPPA, THETA, XI, RHO, LAM, MU_J, SIG_J,
                    OptionType.CALL), abs=1e-9)


def test_bad_lambda_raises():
    with pytest.raises(ValueError):
        bates_greeks(S, K, T, R, V0, KAPPA, THETA, XI, RHO, -1.0, MU_J, SIG_J)
