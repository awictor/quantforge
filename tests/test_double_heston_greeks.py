"""Greeks of a double-Heston option (double_heston_greeks)."""

import pytest

from quantforge import double_heston_greeks, double_heston_price, OptionType


S, K, T, R = 100.0, 100.0, 1.0, 0.03
P = dict(v01=0.03, kappa1=1.5, theta1=0.03, xi1=0.4, rho1=-0.6,
         v02=0.02, kappa2=0.5, theta2=0.02, xi2=0.3, rho2=-0.5)


def test_delta_matches_finite_difference():
    g = double_heston_greeks(S, K, T, R, **P, option_type=OptionType.CALL)
    h = 0.5
    fd = (double_heston_price(S + h, K, T, R, **P)
          - double_heston_price(S - h, K, T, R, **P)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-3)


def test_call_greek_signs():
    g = double_heston_greeks(S, K, T, R, **P, option_type=OptionType.CALL)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0


def test_both_variance_vegas_positive():
    g = double_heston_greeks(S, K, T, R, **P, option_type=OptionType.CALL)
    assert g["vega_v01"] > 0.0
    assert g["vega_v02"] > 0.0


def test_price_field_matches_price():
    g = double_heston_greeks(S, K, T, R, **P, option_type=OptionType.CALL)
    assert g["price"] == pytest.approx(
        double_heston_price(S, K, T, R, **P, option_type=OptionType.CALL),
        abs=1e-9)


def test_bad_variance_raises():
    bad = dict(P)
    bad["v01"] = -0.01
    with pytest.raises(ValueError):
        double_heston_greeks(S, K, T, R, **bad)
