"""Levy-model calibration to a market smile: synthetic-parameter recovery."""

import pytest

from quantforge import (
    calibrate_levy_smile,
    variance_gamma_smile,
    nig_smile,
    meixner_smile,
    cgmy_smile,
)


S, T, R = 100.0, 0.5, 0.03
STRIKES = [80, 90, 100, 110, 120]


def _vols(smile):
    return [iv for _, iv in smile]


@pytest.mark.slow
def test_recovers_vg_parameters():
    true = (0.2, 0.35, -0.25)
    mkt = _vols(variance_gamma_smile(S, STRIKES, T, R, *true))
    p, rmse = calibrate_levy_smile("vg", S, T, R, STRIKES, mkt)
    assert rmse < 1e-3
    assert p[0] == pytest.approx(true[0], abs=1e-2)
    assert p[2] == pytest.approx(true[2], abs=2e-2)


@pytest.mark.slow
def test_recovers_nig_parameters():
    true = (18.0, -6.0, 0.55)
    mkt = _vols(nig_smile(S, STRIKES, T, R, *true))
    p, rmse = calibrate_levy_smile("nig", S, T, R, STRIKES, mkt)
    assert rmse < 1e-3
    assert p[1] == pytest.approx(true[1], abs=0.2)  # asymmetry


@pytest.mark.slow
def test_recovers_meixner_parameters():
    true = (0.35, -0.4, 0.6)
    mkt = _vols(meixner_smile(S, STRIKES, T, R, *true))
    p, rmse = calibrate_levy_smile("meixner", S, T, R, STRIKES, mkt)
    assert rmse < 1e-3


@pytest.mark.slow
def test_recovers_cgmy_fit_quality():
    # CGMY has four parameters and is only weakly identified from five strikes;
    # require a good *fit* (low RMSE) rather than exact parameter recovery.
    true = (4.0, 5.0, 10.0, 0.6)
    mkt = _vols(cgmy_smile(S, STRIKES, T, R, *true))
    _p, rmse = calibrate_levy_smile("cgmy", S, T, R, STRIKES, mkt, max_iter=8000)
    assert rmse < 5e-3


@pytest.mark.slow
def test_cross_model_fit_is_close():
    # A NIG fit to a VG-generated smile should be close but not exact.
    mkt = _vols(variance_gamma_smile(S, STRIKES, T, R, 0.2, 0.4, -0.3))
    _p, rmse = calibrate_levy_smile("nig", S, T, R, STRIKES, mkt)
    assert rmse < 1e-2


def test_unknown_model_raises():
    with pytest.raises(ValueError):
        calibrate_levy_smile("heston", S, T, R, STRIKES, [0.2] * 5)


def test_too_few_points_raises():
    with pytest.raises(ValueError):
        calibrate_levy_smile("vg", S, T, R, [100, 110], [0.2, 0.19])
