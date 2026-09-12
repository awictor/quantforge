"""Forecast-accuracy metrics."""

import random

import pytest

from quantforge import mae, rmse, mape, smape, mase


def test_perfect_forecast_zero_error():
    a = [10, 12, 14, 13, 15]
    assert mae(a, a) == 0.0
    assert rmse(a, a) == 0.0
    assert mape(a, a) == 0.0
    assert smape(a, a) == 0.0


def test_known_values():
    a, f = [10, 20], [12, 18]
    assert mae(a, f) == 2.0
    assert rmse(a, f) == 2.0
    assert abs(mape(a, f) - 0.15) < 1e-12


def test_rmse_at_least_mae():
    rng = random.Random(1)
    a = [rng.gauss(0, 1) for _ in range(100)]
    f = [rng.gauss(0, 1) for _ in range(100)]
    assert rmse(a, f) >= mae(a, f)


def test_mape_scale_invariant():
    a, f = [10, 20], [12, 18]
    scaled_a = [x * 100 for x in a]
    scaled_f = [x * 100 for x in f]
    assert abs(mape(a, f) - mape(scaled_a, scaled_f)) < 1e-12


def test_smape_bounded():
    rng = random.Random(1)
    a = [rng.gauss(0, 1) for _ in range(100)]
    f = [rng.gauss(0, 1) for _ in range(100)]
    assert 0.0 <= smape(a, f) <= 2.0


def test_mase_naive_is_one():
    train = [10, 11, 12, 13, 14, 15]
    act = [16, 17, 18]
    naive_fc = [15, 16, 17]                 # each off by 1, as is the train diff
    assert abs(mase(act, naive_fc, train, 1) - 1.0) < 1e-9


def test_mase_below_one_for_perfect_forecast():
    train = [10, 11, 12, 13, 14, 15]
    act = [16, 17, 18]
    assert mase(act, act, train, 1) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        mae([1], [1, 2])
    with pytest.raises(ValueError):
        mape([0, 1], [1, 1])                 # zero actual
    with pytest.raises(ValueError):
        mase([1, 2], [1, 2], [5], 1)         # train shorter than lag
