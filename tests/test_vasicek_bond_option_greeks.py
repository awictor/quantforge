"""Greeks of a Vasicek zero-coupon-bond option (vasicek_bond_option_greeks)."""

import pytest

from quantforge import vasicek_bond_option_greeks, OptionType
from quantforge.vasicek import bond_option


R0, TO, TB, K = 0.03, 1.0, 5.0, 0.85
KAPPA, THETA, SIG = 0.3, 0.04, 0.01


def test_rho_r_matches_finite_difference():
    g = vasicek_bond_option_greeks(R0, TO, TB, K, KAPPA, THETA, SIG,
                                   OptionType.CALL)
    h = 1e-5
    fd = (bond_option(R0 + h, TO, TB, K, KAPPA, THETA, SIG, OptionType.CALL)
          - bond_option(R0 - h, TO, TB, K, KAPPA, THETA, SIG,
                        OptionType.CALL)) / (2 * h)
    assert g["rho_r"] == pytest.approx(fd, abs=1e-3)


def test_call_rho_r_negative():
    # A bond call loses value as the short rate rises.
    g = vasicek_bond_option_greeks(R0, TO, TB, K, KAPPA, THETA, SIG,
                                   OptionType.CALL)
    assert g["rho_r"] < 0.0


def test_put_rho_r_positive():
    g = vasicek_bond_option_greeks(R0, TO, TB, K, KAPPA, THETA, SIG,
                                   OptionType.PUT)
    assert g["rho_r"] > 0.0


def test_vega_positive():
    g = vasicek_bond_option_greeks(R0, TO, TB, K, KAPPA, THETA, SIG,
                                   OptionType.CALL)
    assert g["vega"] > 0.0


def test_price_field_matches_bond_option():
    g = vasicek_bond_option_greeks(R0, TO, TB, K, KAPPA, THETA, SIG,
                                   OptionType.CALL)
    assert g["price"] == pytest.approx(
        bond_option(R0, TO, TB, K, KAPPA, THETA, SIG, OptionType.CALL),
        abs=1e-12)


def test_bad_times_raise():
    with pytest.raises(ValueError):
        vasicek_bond_option_greeks(R0, 5.0, 1.0, K, KAPPA, THETA, SIG)
