"""Tests for the SVI volatility smile and its calibration."""

import math

import pytest

from quantforge import SVIParams, calibrate_svi
from quantforge.optimize import nelder_mead


# --- Nelder-Mead sanity ---
def test_nelder_mead_quadratic():
    # Minimum of (x-3)^2 + (y+1)^2 is at (3, -1), value 0.
    f = lambda p: (p[0] - 3) ** 2 + (p[1] + 1) ** 2
    x, fx = nelder_mead(f, [0.0, 0.0], tol=1e-14)
    assert x[0] == pytest.approx(3.0, abs=1e-5)
    assert x[1] == pytest.approx(-1.0, abs=1e-5)
    assert fx == pytest.approx(0.0, abs=1e-9)


def test_nelder_mead_rosenbrock():
    # Rosenbrock minimum at (1, 1).
    f = lambda p: (1 - p[0]) ** 2 + 100 * (p[1] - p[0] ** 2) ** 2
    x, fx = nelder_mead(f, [-1.2, 1.0], max_iter=5000, tol=1e-16)
    assert x[0] == pytest.approx(1.0, abs=1e-3)
    assert x[1] == pytest.approx(1.0, abs=1e-3)


# --- SVI formula ---
def test_svi_total_variance_at_min():
    # At k = m the sqrt term is s, so w(m) = a + b*s (rho term vanishes).
    p = SVIParams(a=0.04, b=0.1, rho=-0.3, m=0.0, s=0.1)
    assert p.total_variance(0.0) == pytest.approx(0.04 + 0.1 * 0.1)


def test_svi_implied_vol_positive():
    p = SVIParams(a=0.04, b=0.4, rho=-0.4, m=0.0, s=0.2)
    for k in (-0.5, -0.1, 0.0, 0.1, 0.5):
        assert p.implied_vol(k, t=1.0) > 0


def test_svi_wing_arbitrage_flag():
    ok = SVIParams(a=0.04, b=0.5, rho=-0.3, m=0.0, s=0.1)
    bad = SVIParams(a=0.04, b=2.0, rho=-0.5, m=0.0, s=0.1)  # b*(1+|rho|)=3 > 2
    assert ok.is_arbitrage_free_wings()
    assert not bad.is_arbitrage_free_wings()


# --- Calibration recovers a known surface ---
def test_calibrate_recovers_params():
    true = SVIParams(a=0.03, b=0.25, rho=-0.4, m=0.02, s=0.15)
    ks = [-0.4, -0.25, -0.1, 0.0, 0.1, 0.25, 0.4]
    tv = [true.total_variance(k) for k in ks]
    fitted, rmse = calibrate_svi(ks, tv)
    assert rmse < 1e-5
    # Fitted surface should reproduce total variance at every quote.
    for k in ks:
        assert fitted.total_variance(k) == pytest.approx(true.total_variance(k), abs=1e-4)


def test_calibrate_fits_noisy_smile():
    # A typical downward-skewed equity smile in vol space.
    t = 0.5
    ks = [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]
    vols = [0.28, 0.24, 0.21, 0.20, 0.205, 0.22, 0.245]
    tv = [(v * v) * t for v in vols]
    fitted, rmse = calibrate_svi(ks, tv)
    assert rmse < 5e-3
    # Recovered vols should be close to the inputs.
    for k, v in zip(ks, vols):
        assert fitted.implied_vol(k, t) == pytest.approx(v, abs=0.02)


def test_calibrate_requires_five_quotes():
    with pytest.raises(ValueError):
        calibrate_svi([0.0, 0.1], [0.04, 0.05])
