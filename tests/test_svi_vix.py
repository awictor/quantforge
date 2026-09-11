"""VIX / SVIX indices implied by an SVI slice (svi_vix, svi_svix)."""

import pytest

from quantforge import svi_vix, svi_svix, SVIParams


S0, T, R = 100.0, 1.0, 0.03


def _flat(sig):
    return SVIParams(a=sig * sig * T, b=0.0, rho=0.0, m=0.0, s=0.1)


def test_flat_vix_is_hundred_sigma():
    v = svi_vix(_flat(0.2), S0, T, R, n_strikes=401, width=8.0)
    assert v == pytest.approx(20.0, abs=0.2)


def test_flat_svix_is_hundred_sigma():
    v = svi_svix(_flat(0.2), S0, T, R, n_strikes=401, width=8.0)
    assert v == pytest.approx(20.0, abs=0.3)


def test_skew_lifts_vix_above_atm():
    p = SVIParams(a=0.03, b=0.1, rho=-0.5, m=0.0, s=0.2)
    atm = 100.0 * p.implied_vol(0.0, T)
    assert svi_vix(p, S0, T, R) > atm


def test_both_positive():
    p = SVIParams(a=0.03, b=0.1, rho=-0.5, m=0.0, s=0.2)
    assert svi_vix(p, S0, T, R) > 0.0
    assert svi_svix(p, S0, T, R) > 0.0


def test_higher_convexity_raises_vix():
    lo = SVIParams(a=0.04, b=0.05, rho=0.0, m=0.0, s=0.3)
    hi = SVIParams(a=0.04, b=0.15, rho=0.0, m=0.0, s=0.3)
    assert svi_vix(hi, S0, T, R) > svi_vix(lo, S0, T, R)
