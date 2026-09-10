"""Greeks of a NIG (normal inverse Gaussian) option (nig_greeks)."""

import pytest

from quantforge import nig_greeks, nig_price, OptionType


S, K, T, R = 100.0, 100.0, 1.0, 0.05
ALPHA, BETA, DELTA = 15.0, -5.0, 10.0


def test_delta_matches_finite_difference():
    g = nig_greeks(S, K, T, R, ALPHA, BETA, DELTA, OptionType.CALL)
    h = 0.01
    fd = (nig_price(S + h, K, T, R, ALPHA, BETA, DELTA)
          - nig_price(S - h, K, T, R, ALPHA, BETA, DELTA)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-4)


def test_call_greek_signs():
    g = nig_greeks(S, K, T, R, ALPHA, BETA, DELTA, OptionType.CALL)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0


def test_put_delta_negative():
    g = nig_greeks(S, K, T, R, ALPHA, BETA, DELTA, OptionType.PUT)
    assert g["delta"] < 0.0


def test_parameter_sensitivities_finite():
    g = nig_greeks(S, K, T, R, ALPHA, BETA, DELTA, OptionType.CALL)
    assert g["d_alpha"] == g["d_alpha"]  # not NaN
    assert g["d_beta"] == g["d_beta"]


def test_price_field_matches_price():
    g = nig_greeks(S, K, T, R, ALPHA, BETA, DELTA, OptionType.CALL)
    assert g["price"] == pytest.approx(
        nig_price(S, K, T, R, ALPHA, BETA, DELTA, OptionType.CALL), abs=1e-9)


def test_bad_beta_raises():
    with pytest.raises(ValueError):
        nig_greeks(S, K, T, R, ALPHA, 20.0, DELTA)  # |beta| >= alpha
