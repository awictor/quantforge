"""Greeks of a quanto option (quanto_option_greeks)."""

import pytest

from quantforge import (
    quanto_option_greeks,
    quanto_option,
    delta as bsm_delta,
    OptionType,
)


S, K, T = 100.0, 100.0, 1.0
RD, RF = 0.04, 0.06
SA, SFX, RHO, Q = 0.2, 0.1, 0.3, 0.0


def test_delta_matches_finite_difference():
    g = quanto_option_greeks(S, K, T, RD, RF, SA, SFX, RHO, Q, OptionType.CALL)
    h = 0.01
    fd = (quanto_option(S + h, K, T, RD, RF, SA, SFX, RHO, Q)
          - quanto_option(S - h, K, T, RD, RF, SA, SFX, RHO, Q)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-5)


def test_zero_correlation_delta_matches_plain_bsm():
    # rho = 0 removes the quanto adjustment: carry is b = r_foreign - q.
    g = quanto_option_greeks(S, K, T, RD, RF, SA, SFX, 0.0, Q, OptionType.CALL)
    bd = bsm_delta(S, K, T, RD, SA, OptionType.CALL, b=RF - Q)
    assert g["delta"] == pytest.approx(bd, abs=1e-9)


def test_correlation_vega_negative_for_positive_rho_call():
    # Higher rho lowers the quanto carry b_q, so the call is worth less.
    g = quanto_option_greeks(S, K, T, RD, RF, SA, SFX, RHO, Q, OptionType.CALL)
    assert g["corr_vega"] < 0.0


def test_call_gamma_and_asset_vega_positive():
    g = quanto_option_greeks(S, K, T, RD, RF, SA, SFX, RHO, Q, OptionType.CALL)
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0


def test_price_field_matches_quanto_option():
    g = quanto_option_greeks(S, K, T, RD, RF, SA, SFX, RHO, Q, OptionType.CALL)
    assert g["price"] == pytest.approx(
        quanto_option(S, K, T, RD, RF, SA, SFX, RHO, Q, OptionType.CALL),
        abs=1e-12)


def test_bad_rho_raises():
    with pytest.raises(ValueError):
        quanto_option_greeks(S, K, T, RD, RF, SA, SFX, 1.5, Q)
