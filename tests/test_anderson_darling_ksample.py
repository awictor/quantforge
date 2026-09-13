"""Anderson-Darling k-sample test (Scholz-Stephens)."""

import random

import pytest

from quantforge import anderson_darling_ksample


def test_scholz_stephens_worked_example():
    s1 = [38.7, 41.5, 43.8, 44.5, 45.5, 46.0, 47.7, 58.0]
    s2 = [39.2, 39.3, 39.7, 41.4, 41.8, 42.9, 43.3, 45.8]
    s3 = [34.0, 35.0, 39.0, 40.0, 43.0, 43.0, 44.0, 45.0]
    s4 = [34.0, 34.8, 34.8, 35.4, 37.2, 37.8, 41.2, 42.8]
    r = anderson_darling_ksample(s1, s2, s3, s4)
    assert abs(r["a2k"] - 8.36) < 0.1              # matches the published statistic
    assert r["p_value"] < 0.01


def test_same_distribution_not_significant():
    rng = random.Random(3)
    groups = [[rng.gauss(0, 1) for _ in range(40)] for _ in range(3)]
    assert anderson_darling_ksample(*groups)["p_value"] > 0.05


def test_different_means_significant():
    rng = random.Random(3)
    groups = [[rng.gauss(m, 1) for _ in range(40)] for m in (0, 2, 4)]
    assert anderson_darling_ksample(*groups)["p_value"] < 0.05


def test_null_calibration():
    rng = random.Random(5)
    rej = 0
    N = 500
    for _ in range(N):
        g = [[rng.gauss(0, 1) for _ in range(25)] for _ in range(3)]
        if anderson_darling_ksample(*g)["p_value"] < 0.05:
            rej += 1
    assert 0.02 <= rej / N <= 0.10


def test_two_sample():
    r = anderson_darling_ksample([1, 2, 3, 4, 5], [6, 7, 8, 9, 10])
    assert r["p_value"] < 0.05


def test_validation():
    with pytest.raises(ValueError):
        anderson_darling_ksample([1, 2, 3])
    with pytest.raises(ValueError):
        anderson_darling_ksample([5, 5], [5, 5])       # all identical
