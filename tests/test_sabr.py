"""Tests for the SABR model: Hagan vol formula and calibration."""

import math

import pytest

from quantforge import SABRParams, sabr_vol, calibrate_sabr


BASE = dict(F=100.0, t=1.0, alpha=0.2, beta=0.5, rho=-0.3, nu=0.4)


def test_atm_limit_is_continuous():
    # The ATM branch (F==K) must match the K->F limit of the general branch.
    F = BASE["F"]
    atm = sabr_vol(F, F, BASE["t"], BASE["alpha"], BASE["beta"], BASE["rho"], BASE["nu"])
    near = sabr_vol(F, F + 1e-6, BASE["t"], BASE["alpha"], BASE["beta"],
                    BASE["rho"], BASE["nu"])
    assert atm == pytest.approx(near, abs=1e-6)


def test_vol_positive_across_strikes():
    for K in (60, 80, 100, 120, 150):
        v = sabr_vol(BASE["F"], K, BASE["t"], BASE["alpha"], BASE["beta"],
                     BASE["rho"], BASE["nu"])
        assert v > 0


def test_zero_volvol_is_flat_cev_level():
    # nu = 0 removes vol-of-vol; the ATM vol is alpha / F^{1-beta} plus the
    # small beta-backbone time correction that SABR retains even at nu=0.
    F, alpha, beta = 100.0, 0.2, 0.5
    v = sabr_vol(F, F, 1.0, alpha, beta, rho=0.0, nu=0.0)
    base = alpha / F ** (1 - beta)
    backbone = (1 - beta) ** 2 / 24.0 * alpha * alpha / (F ** (1 - beta)) ** 2
    assert v == pytest.approx(base * (1.0 + backbone * 1.0), abs=1e-12)


def test_negative_rho_gives_downward_skew():
    # Equity-style negative correlation -> higher vol for low strikes.
    lo = sabr_vol(100, 80, 1.0, 0.2, 0.5, rho=-0.5, nu=0.5)
    hi = sabr_vol(100, 120, 1.0, 0.2, 0.5, rho=-0.5, nu=0.5)
    assert lo > hi


def test_calibration_recovers_known_smile():
    true = SABRParams(alpha=0.25, beta=0.5, rho=-0.4, nu=0.6)
    F, t = 100.0, 1.0
    strikes = [70, 85, 100, 115, 130]
    vols = [sabr_vol(F, K, t, true.alpha, true.beta, true.rho, true.nu) for K in strikes]

    fitted, rmse = calibrate_sabr(F, t, strikes, vols, beta=0.5)
    assert rmse < 1e-5
    # Recovered smile reproduces every market vol.
    for K, v in zip(strikes, vols):
        model = sabr_vol(F, K, t, fitted.alpha, fitted.beta, fitted.rho, fitted.nu)
        assert model == pytest.approx(v, abs=1e-4)


def test_calibration_fits_noisy_smile():
    F, t = 100.0, 0.5
    strikes = [80, 90, 100, 110, 120]
    market = [0.26, 0.235, 0.22, 0.225, 0.24]  # a realistic smile
    fitted, rmse = calibrate_sabr(F, t, strikes, market, beta=0.5)
    assert rmse < 5e-3
    for K, v in zip(strikes, market):
        model = sabr_vol(F, K, t, fitted.alpha, fitted.beta, fitted.rho, fitted.nu)
        assert model == pytest.approx(v, abs=0.02)


def test_calibration_needs_three_points():
    with pytest.raises(ValueError):
        calibrate_sabr(100, 1.0, [90, 110], [0.2, 0.21])


def test_invalid_inputs_rejected():
    with pytest.raises(ValueError):
        sabr_vol(-1, 100, 1.0, 0.2, 0.5, -0.3, 0.4)
    with pytest.raises(ValueError):
        sabr_vol(100, 100, 1.0, 0.0, 0.5, -0.3, 0.4)  # alpha must be > 0


def test_zero_vol_of_vol_no_divide_by_zero():
    # nu = 0 makes z = 0; the z/x(z) ratio must take its limit of 1 rather than
    # divide by zero. The result matches the nu -> 0 limit.
    off_atm = sabr_vol(100, 110, 1.0, 0.2, 0.5, -0.3, 0.0)
    limit = sabr_vol(100, 110, 1.0, 0.2, 0.5, -0.3, 1e-8)
    assert off_atm == pytest.approx(limit, abs=1e-7)


def test_zero_vol_of_vol_still_has_beta_skew():
    # With nu = 0 the smile is flat in vol-of-vol but still skewed by beta < 1.
    lo = sabr_vol(100, 80, 1.0, 0.2, 0.5, -0.3, 0.0)
    hi = sabr_vol(100, 120, 1.0, 0.2, 0.5, -0.3, 0.0)
    assert lo > hi  # downward skew from beta


def test_zero_vol_of_vol_atm_is_alpha_over_fbeta():
    # ATM with nu = 0, rho = 0: sigma ~ alpha / F^{1-beta} plus a tiny
    # beta^2 alpha^2 time correction, so it sits just above the leading term.
    atm = sabr_vol(100, 100, 1.0, 0.2, 0.5, 0.0, 0.0)
    assert atm == pytest.approx(0.2 / 100 ** 0.5, abs=1e-4)
