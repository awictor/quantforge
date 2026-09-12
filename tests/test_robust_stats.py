"""Robust scale and location estimators."""

import random
import statistics

import pytest

from quantforge import (
    median_absolute_deviation, interquartile_range, winsorize, trimmed_mean,
)


def test_mad_consistent_with_std_on_gaussian():
    rng = random.Random(1)
    x = [rng.gauss(0, 2.0) for _ in range(5000)]
    assert abs(median_absolute_deviation(x) - 2.0) < 0.15


def test_iqr_scaled_consistent_with_std():
    rng = random.Random(1)
    x = [rng.gauss(0, 2.0) for _ in range(5000)]
    assert abs(interquartile_range(x, scale=True) - 2.0) < 0.15


def test_constant_data_zero_scale():
    assert median_absolute_deviation([5, 5, 5, 5]) == 0.0
    assert interquartile_range([5, 5, 5, 5]) == 0.0


def test_mad_robust_where_std_explodes():
    rng = random.Random(1)
    x = [rng.gauss(0, 2.0) for _ in range(100)] + [1e6] * 10
    assert median_absolute_deviation(x) < 5.0
    assert statistics.pstdev(x) > 1000.0


def test_winsorize_clips_and_preserves_length():
    w = winsorize([1, 2, 3, 4, 5, 6, 7, 8, 9, 100], limit=0.1)
    assert len(w) == 10
    assert max(w) < 100


def test_trimmed_mean_robust():
    data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 1000]
    assert abs(trimmed_mean(data, 0.1) - 5.5) < 1.0
    assert abs(trimmed_mean(data, 0.0) - sum(data) / 10) < 1e-12


def test_unscaled_mad_is_raw_median_deviation():
    # median of |x - median| for a symmetric set is exact.
    x = [1, 2, 3, 4, 5]        # median 3; deviations 2,1,0,1,2 -> median 1
    assert median_absolute_deviation(x, scale=False) == 1.0


def test_validation():
    with pytest.raises(ValueError):
        median_absolute_deviation([])
    with pytest.raises(ValueError):
        interquartile_range([1.0])
    with pytest.raises(ValueError):
        winsorize([1, 2, 3], limit=0.5)
    with pytest.raises(ValueError):
        trimmed_mean([1, 2, 3], proportion=0.6)
