"""Equal-variance tests: Levene / Brown-Forsythe and Bartlett."""

import random

import pytest

from quantforge import levene_test, bartlett_test


def test_equal_variance_not_significant():
    rng = random.Random(3)
    groups = [[rng.gauss(0, 1) for _ in range(50)] for _ in range(3)]
    assert levene_test(*groups)["p_value"] > 0.05


def test_unequal_variance_significant():
    rng = random.Random(3)
    groups = [[rng.gauss(0, s) for _ in range(80)] for s in (1, 3, 6)]
    assert levene_test(*groups)["p_value"] < 0.01
    assert bartlett_test(*groups)["p_value"] < 0.01


def test_levene_mean_and_median_variants():
    rng = random.Random(9)
    groups = [[rng.gauss(0, 1) for _ in range(40)] for _ in range(3)]
    assert levene_test(*groups, center="mean")["p_value"] > 0.05
    assert levene_test(*groups, center="median")["p_value"] > 0.05


def test_levene_null_calibration():
    rng = random.Random(5)
    rej = 0
    N = 500
    for _ in range(N):
        g = [[rng.gauss(0, 1) for _ in range(25)] for _ in range(3)]
        if levene_test(*g)["p_value"] < 0.05:
            rej += 1
    assert 0.02 <= rej / N <= 0.08


def test_bartlett_two_group():
    r = bartlett_test([1, 2, 3, 4, 5], [1, 3, 5, 7, 9])
    assert r["df"] == 1
    assert 0.0 <= r["p_value"] <= 1.0


def test_validation():
    with pytest.raises(ValueError):
        levene_test([1, 2, 3])
    with pytest.raises(ValueError):
        levene_test([1, 2], [3, 4], center="bad")
    with pytest.raises(ValueError):
        bartlett_test([1], [2, 3])                # a group with < 2 obs
