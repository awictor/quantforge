"""Rate sensitivities of Vasicek and Ho-Lee zero-coupon bonds."""

import pytest

from quantforge import (
    vasicek_bond_greeks, holee_bond_greeks,
)
from quantforge.vasicek import zero_coupon_bond as vasicek_zcb
from quantforge.holee import holee_zero_coupon_bond


VK, VT, KAPPA, VTHETA, VSIG = 0.03, 5.0, 0.3, 0.04, 0.02
HR, HT, HTHETA, HSIG = 0.03, 5.0, 0.01, 0.02


def test_vasicek_rho_matches_finite_difference():
    g = vasicek_bond_greeks(VK, VT, KAPPA, VTHETA, VSIG)
    h = 1e-6
    fd = (vasicek_zcb(VK + h, VT, KAPPA, VTHETA, VSIG)
          - vasicek_zcb(VK - h, VT, KAPPA, VTHETA, VSIG)) / (2 * h)
    assert g["rho_r"] == pytest.approx(fd, abs=1e-6)


def test_vasicek_signs_and_convexity():
    g = vasicek_bond_greeks(VK, VT, KAPPA, VTHETA, VSIG)
    assert g["rho_r"] < 0.0
    assert g["gamma_r"] > 0.0
    assert g["convexity"] == pytest.approx(g["duration"] ** 2, abs=1e-9)


def test_holee_rho_matches_finite_difference():
    g = holee_bond_greeks(HR, HT, HTHETA, HSIG)
    h = 1e-6
    fd = (holee_zero_coupon_bond(HR + h, HT, HTHETA, HSIG)
          - holee_zero_coupon_bond(HR - h, HT, HTHETA, HSIG)) / (2 * h)
    assert g["rho_r"] == pytest.approx(fd, abs=1e-6)


def test_holee_duration_equals_maturity():
    # A Ho-Lee zero-coupon bond has rate duration exactly equal to its maturity.
    g = holee_bond_greeks(HR, HT, HTHETA, HSIG)
    assert g["duration"] == pytest.approx(HT, abs=1e-12)
    assert g["convexity"] == pytest.approx(HT * HT, abs=1e-12)


def test_zero_maturity_flat():
    gv = vasicek_bond_greeks(VK, 0.0, KAPPA, VTHETA, VSIG)
    assert gv["price"] == 1.0 and gv["duration"] == 0.0
    gh = holee_bond_greeks(HR, 0.0, HTHETA, HSIG)
    assert gh["price"] == 1.0 and gh["duration"] == 0.0


def test_bad_maturity_raises():
    with pytest.raises(ValueError):
        vasicek_bond_greeks(VK, -1.0, KAPPA, VTHETA, VSIG)
