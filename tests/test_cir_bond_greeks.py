"""Rate sensitivities of a CIR zero-coupon bond (cir_bond_greeks)."""

import pytest

from quantforge import cir_bond_greeks, cir_zero_coupon_bond


R0, T, KAPPA, THETA, SIG = 0.03, 5.0, 0.3, 0.04, 0.02


def test_rho_r_matches_finite_difference():
    g = cir_bond_greeks(R0, T, KAPPA, THETA, SIG)
    h = 1e-6
    fd = (cir_zero_coupon_bond(R0 + h, T, KAPPA, THETA, SIG)
          - cir_zero_coupon_bond(R0 - h, T, KAPPA, THETA, SIG)) / (2 * h)
    assert g["rho_r"] == pytest.approx(fd, abs=1e-6)


def test_gamma_r_matches_finite_difference():
    g = cir_bond_greeks(R0, T, KAPPA, THETA, SIG)
    h = 1e-4
    fg = (cir_zero_coupon_bond(R0 + h, T, KAPPA, THETA, SIG)
          - 2 * cir_zero_coupon_bond(R0, T, KAPPA, THETA, SIG)
          + cir_zero_coupon_bond(R0 - h, T, KAPPA, THETA, SIG)) / (h * h)
    assert g["gamma_r"] == pytest.approx(fg, rel=1e-3)


def test_signs_and_duration_convexity():
    g = cir_bond_greeks(R0, T, KAPPA, THETA, SIG)
    assert g["rho_r"] < 0.0          # bond falls as the short rate rises
    assert g["gamma_r"] > 0.0        # convex in the rate
    assert g["duration"] > 0.0
    assert g["convexity"] == pytest.approx(g["duration"] ** 2, abs=1e-9)


def test_duration_relates_to_rho():
    g = cir_bond_greeks(R0, T, KAPPA, THETA, SIG)
    assert g["duration"] == pytest.approx(-g["rho_r"] / g["price"], abs=1e-9)


def test_zero_maturity_is_flat():
    g = cir_bond_greeks(R0, 0.0, KAPPA, THETA, SIG)
    assert g["price"] == 1.0
    assert g["rho_r"] == 0.0
    assert g["duration"] == 0.0


def test_bad_params_raise():
    with pytest.raises(ValueError):
        cir_bond_greeks(-0.01, T, KAPPA, THETA, SIG)
