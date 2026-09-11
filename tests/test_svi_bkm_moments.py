"""Risk-neutral moments (BKM) implied by an SVI slice (svi_bkm_moments)."""

import pytest

from quantforge import svi_bkm_moments, SVIParams


S0, T, R = 100.0, 1.0, 0.03


def test_flat_slice_near_symmetric():
    sig = 0.2
    flat = SVIParams(a=sig * sig * T, b=0.0, rho=0.0, m=0.0, s=0.1)
    var, skew, kurt = svi_bkm_moments(flat, S0, T, R)
    assert var == pytest.approx(sig * sig, abs=2e-3)
    assert abs(skew) < 0.05
    assert abs(kurt) < 0.1


def test_negative_rho_gives_negative_skew():
    p = SVIParams(a=0.03, b=0.1, rho=-0.6, m=0.0, s=0.2)
    _var, skew, _kurt = svi_bkm_moments(p, S0, T, R)
    assert skew < 0.0


def test_positive_rho_gives_positive_skew():
    p = SVIParams(a=0.03, b=0.1, rho=0.6, m=0.0, s=0.2)
    _var, skew, _kurt = svi_bkm_moments(p, S0, T, R)
    assert skew > 0.0


def test_convex_slice_positive_excess_kurtosis():
    p = SVIParams(a=0.03, b=0.15, rho=0.0, m=0.0, s=0.15)
    _var, _skew, kurt = svi_bkm_moments(p, S0, T, R)
    assert kurt > 0.0


def test_variance_positive():
    p = SVIParams(a=0.03, b=0.1, rho=-0.3, m=0.0, s=0.2)
    var, _skew, _kurt = svi_bkm_moments(p, S0, T, R)
    assert var > 0.0
