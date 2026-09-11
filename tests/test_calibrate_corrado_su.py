"""Calibration of Corrado-Su (sigma, skew, excess kurtosis) to call prices."""

import pytest

from quantforge import corrado_su_call, calibrate_corrado_su
from quantforge.bsm import call_price


S, T, R = 100.0, 0.5, 0.03
KS = [80.0, 90.0, 95.0, 100.0, 105.0, 110.0, 120.0]


def test_recovers_known_parameters():
    sigma, skew, kurt = 0.22, -0.6, 1.5
    prices = [corrado_su_call(S, K, T, R, sigma, skew, kurt) for K in KS]
    fs, fsk, fku, rmse = calibrate_corrado_su(S, T, R, KS, prices)
    assert fs == pytest.approx(sigma, abs=1e-3)
    assert fsk == pytest.approx(skew, abs=1e-2)
    assert fku == pytest.approx(kurt, abs=1e-2)
    assert rmse < 1e-4


def test_flat_bs_calibrates_to_zero_moments():
    prices = [call_price(S, K, T, R, 0.25) for K in KS]
    fs, fsk, fku, rmse = calibrate_corrado_su(S, T, R, KS, prices)
    assert fs == pytest.approx(0.25, abs=1e-3)
    assert abs(fsk) < 1e-2
    assert abs(fku) < 1e-2
    assert rmse < 1e-4


def test_reprices_within_rmse():
    sigma, skew, kurt = 0.3, 0.4, 2.0
    prices = [corrado_su_call(S, K, T, R, sigma, skew, kurt) for K in KS]
    fs, fsk, fku, rmse = calibrate_corrado_su(S, T, R, KS, prices)
    for K, mkt in zip(KS, prices):
        model = corrado_su_call(S, K, T, R, fs, fsk, fku)
        assert model == pytest.approx(mkt, abs=5e-3)


def test_too_few_quotes_raises():
    with pytest.raises(ValueError):
        calibrate_corrado_su(S, T, R, [100.0, 105.0], [3.0, 1.5])
