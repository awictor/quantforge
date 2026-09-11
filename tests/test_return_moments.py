"""Sample skewness, kurtosis, and Jarque-Bera (perfmetrics module)."""

import random

import pytest

from quantforge import sample_skewness, sample_kurtosis, jarque_bera


def test_normal_sample_near_zero():
    rng = random.Random(3)
    norm = [rng.gauss(0.0, 1.0) for _ in range(20000)]
    assert abs(sample_skewness(norm)) < 0.1
    assert abs(sample_kurtosis(norm)) < 0.15
    assert jarque_bera(norm) < 6.0   # below the 5% chi-sq(2) critical value


def test_symmetric_zero_skew():
    assert sample_skewness([-2, -1, 0, 1, 2]) == pytest.approx(0.0, abs=1e-12)


def test_right_skew_positive():
    assert sample_skewness([0, 0, 0, 0, 10]) > 0.0


def test_fat_tails_positive_excess_kurtosis():
    fat = [0.0] * 100 + [10.0, -10.0]
    assert sample_kurtosis(fat) > 0.0


def test_excess_vs_raw_differ_by_three():
    rng = random.Random(5)
    x = [rng.gauss(0.0, 1.0) for _ in range(5000)]
    assert sample_kurtosis(x, excess=False) - sample_kurtosis(x, excess=True) == \
        pytest.approx(3.0, abs=1e-12)


def test_jarque_bera_rejects_skewed():
    skewed = [0.0] * 50 + [20.0]
    assert jarque_bera(skewed) > 6.0


def test_validation():
    with pytest.raises(ValueError):
        sample_skewness([1.0, 1.0, 1.0])   # zero variance
    with pytest.raises(ValueError):
        sample_kurtosis([1.0])             # need >= 2
