"""Corrado-Su implied-vol smile (corrado_su_implied_vol, corrado_su_smile)."""

import pytest

from quantforge import (
    corrado_su_implied_vol, corrado_su_smile, corrado_su_call,
    calibrate_corrado_su,
)


S, T, R = 100.0, 0.5, 0.03
KS = [80.0, 90.0, 100.0, 110.0, 120.0]


def test_flat_when_no_moments():
    for K in KS:
        assert corrado_su_implied_vol(S, K, T, R, 0.25) == pytest.approx(
            0.25, abs=1e-6)


def test_negative_skew_lifts_put_wing():
    _, vols = corrado_su_smile(S, T, R, 0.22, KS, skew=-0.8, excess_kurt=1.0)
    # Low-strike (put) wing richer than the high-strike (call) wing.
    assert vols[0] > vols[-1]
    # Monotone decreasing across strikes for a strong negative skew.
    assert all(vols[i] > vols[i + 1] for i in range(len(vols) - 1))


def test_positive_kurtosis_lifts_both_wings():
    _, vols = corrado_su_smile(S, T, R, 0.22, [80.0, 100.0, 120.0],
                               skew=0.0, excess_kurt=2.0)
    assert vols[0] > vols[1]     # low-strike wing above ATM
    assert vols[2] > vols[1]     # high-strike wing above ATM


def test_implied_vol_reprices_corrado_su():
    sigma, skew, kurt = 0.28, -0.5, 1.8
    K = 95.0
    iv = corrado_su_implied_vol(S, K, T, R, sigma, skew, kurt)
    from quantforge.bsm import call_price
    cs = corrado_su_call(S, K, T, R, sigma, skew, kurt)
    assert call_price(S, K, T, R, iv) == pytest.approx(cs, abs=1e-6)


def test_smile_round_trips_through_calibration():
    # Build a Corrado-Su price surface, calibrate, and confirm the recovered
    # parameters reproduce the same smile.
    sigma, skew, kurt = 0.24, -0.6, 1.2
    prices = [corrado_su_call(S, K, T, R, sigma, skew, kurt) for K in KS]
    fs, fsk, fku, rmse = calibrate_corrado_su(S, T, R, KS, prices)
    _, v_true = corrado_su_smile(S, T, R, sigma, KS, skew=skew, excess_kurt=kurt)
    _, v_fit = corrado_su_smile(S, T, R, fs, KS, skew=fsk, excess_kurt=fku)
    for a, b_ in zip(v_true, v_fit):
        assert a == pytest.approx(b_, abs=1e-3)
