"""Serial-correlation tests: Ljung-Box, Box-Pierce, Durbin-Watson."""

import random

import pytest

from quantforge import ljung_box, box_pierce, durbin_watson
from quantforge.serial_correlation import _chi2_sf


def _white(n, seed):
    rng = random.Random(seed)
    return [rng.gauss(0, 1) for _ in range(n)]


def _ar1(n, phi, seed):
    rng = random.Random(seed)
    y = [0.0]
    for _ in range(n):
        y.append(phi * y[-1] + rng.gauss(0, 1))
    return y[1:]


def test_chi2_survival_reference_points():
    assert abs(_chi2_sf(9.34, 10) - 0.5) < 0.02      # median of chi2(10)
    assert abs(_chi2_sf(3.841, 1) - 0.05) < 0.005    # 95th percentile of chi2(1)


def test_white_noise_fails_to_reject():
    q, p = ljung_box(_white(1000, 1), 10)
    assert p > 0.05


def test_ar1_rejects():
    q, p = ljung_box(_ar1(1000, 0.7, 1), 10)
    assert q > 50.0
    assert p < 1e-6


def test_ljung_box_at_least_box_pierce():
    ar = _ar1(1000, 0.7, 1)
    q_lb, _ = ljung_box(ar, 10)
    q_bp, _ = box_pierce(ar, 10)
    assert q_lb >= q_bp


def test_durbin_watson_regimes():
    assert abs(durbin_watson(_white(1000, 1)) - 2.0) < 0.2
    assert durbin_watson(_ar1(1000, 0.7, 1)) < 1.0        # positive autocorrelation
    assert durbin_watson(_ar1(1000, -0.7, 2)) > 3.0       # negative autocorrelation


def test_validation():
    wn = _white(100, 1)
    with pytest.raises(ValueError):
        ljung_box(wn, 0)
    with pytest.raises(ValueError):
        ljung_box(wn, 100)                # lags >= n
    with pytest.raises(ValueError):
        durbin_watson([1.0])
