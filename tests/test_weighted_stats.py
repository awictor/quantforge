"""Weighted descriptive statistics."""

import random
import statistics

import pytest

from quantforge import (
    weighted_mean,
    weighted_variance,
    weighted_std,
    weighted_quantile,
    weighted_median,
)


def test_equal_weights_match_ordinary():
    rng = random.Random(1)
    for _ in range(500):
        x = [rng.uniform(-10, 10) for _ in range(rng.randint(2, 20))]
        w = [1.0] * len(x)
        assert abs(weighted_mean(x, w) - statistics.mean(x)) < 1e-9
        assert abs(weighted_variance(x, w) - statistics.variance(x)) < 1e-9
        assert abs(weighted_std(x, w) - statistics.stdev(x)) < 1e-9


def test_weighted_mean_known():
    assert weighted_mean([1, 2, 3], [1, 1, 4]) == 2.5


def test_repeat_weight_equivalence():
    rng = random.Random(2)
    for _ in range(200):
        x = [rng.uniform(0, 10) for _ in range(rng.randint(2, 6))]
        w = [rng.randint(1, 4) for _ in x]
        repl = []
        for xi, wi in zip(x, w):
            repl.extend([xi] * wi)
        assert abs(weighted_mean(x, w) - statistics.mean(repl)) < 1e-9
        assert abs(weighted_variance(x, w, unbiased=False) -
                   statistics.pvariance(repl)) < 1e-9


def test_weighted_median():
    assert abs(weighted_median([1, 2, 3, 4], [1, 1, 1, 1]) - 2.5) < 1e-9
    assert weighted_median([1, 2, 3], [1, 1, 10]) > 2.5      # skewed toward 3


def test_weighted_quantile_bounds():
    x = list(range(1, 101))
    w = [1.0] * 100
    assert weighted_quantile(x, w, 0.0) == 1
    assert weighted_quantile(x, w, 1.0) == 100
    assert 49 < weighted_quantile(x, w, 0.5) < 52


def test_validation():
    with pytest.raises(ValueError):
        weighted_mean([1, 2], [1])
    with pytest.raises(ValueError):
        weighted_mean([1, 2], [0, 0])
    with pytest.raises(ValueError):
        weighted_quantile([1, 2], [1, 1], 1.5)
