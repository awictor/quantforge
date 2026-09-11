"""Risk-neutral moments (BKM) implied by a SABR smile (sabr_bkm_moments)."""

import pytest

from quantforge import sabr_bkm_moments


F, T, R = 100.0, 1.0, 0.03


def test_flat_near_symmetric():
    var, skew, kurt = sabr_bkm_moments(F, T, R, 0.2, 1.0, 0.0, 1e-8)
    assert var == pytest.approx(0.04, abs=2e-3)
    assert abs(skew) < 0.05
    assert abs(kurt) < 0.1


def test_negative_rho_negative_skew():
    _v, skew, kurt = sabr_bkm_moments(F, T, R, 0.2, 0.5, -0.6, 0.4)
    assert skew < 0.0
    assert kurt > 0.0


def test_positive_rho_positive_skew():
    _v, skew, _k = sabr_bkm_moments(F, T, R, 0.2, 0.5, 0.6, 0.4)
    assert skew > 0.0


def test_variance_positive():
    var, _s, _k = sabr_bkm_moments(F, T, R, 0.2, 0.5, -0.3, 0.4)
    assert var > 0.0


def test_higher_volvol_raises_kurtosis():
    _v1, _s1, k_lo = sabr_bkm_moments(F, T, R, 0.2, 0.5, 0.0, 0.2)
    _v2, _s2, k_hi = sabr_bkm_moments(F, T, R, 0.2, 0.5, 0.0, 0.6)
    assert k_hi > k_lo
