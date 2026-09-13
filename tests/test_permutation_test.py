"""General permutation tests."""

import random
import statistics

import pytest

from quantforge import permutation_test, paired_permutation_test


def test_detects_mean_difference():
    rng = random.Random(3)
    a = [rng.gauss(0, 1) for _ in range(40)]
    b = [rng.gauss(2, 1) for _ in range(40)]
    assert permutation_test(a, b, n_permutations=2000)["p_value"] < 0.01


def test_same_distribution_not_significant():
    rng = random.Random(3)
    a = [rng.gauss(0, 1) for _ in range(40)]
    b = [rng.gauss(0, 1) for _ in range(40)]
    assert permutation_test(a, b, n_permutations=2000)["p_value"] > 0.05


def test_custom_statistic_median():
    rng = random.Random(5)
    a = [rng.gauss(0, 1) for _ in range(50)]
    b = [rng.gauss(1.5, 1) for _ in range(50)]
    med = lambda a, b: statistics.median(a) - statistics.median(b)
    assert permutation_test(a, b, statistic=med, n_permutations=2000)["p_value"] < 0.01


def test_deterministic():
    rng = random.Random(7)
    a = [rng.gauss(0, 1) for _ in range(30)]
    b = [rng.gauss(0.5, 1) for _ in range(30)]
    r1 = permutation_test(a, b, seed=42, n_permutations=500)
    r2 = permutation_test(a, b, seed=42, n_permutations=500)
    assert r1["p_value"] == r2["p_value"]


def test_null_calibration():
    rng = random.Random(5)
    rej = 0
    N = 500
    for _ in range(N):
        a = [rng.gauss(0, 1) for _ in range(20)]
        b = [rng.gauss(0, 1) for _ in range(20)]
        if permutation_test(a, b, n_permutations=200)["p_value"] < 0.05:
            rej += 1
    assert 0.02 <= rej / N <= 0.08


def test_paired_detects_shift():
    rng = random.Random(7)
    x = [rng.gauss(0, 1) for _ in range(30)]
    y = [xi - 0.8 + rng.gauss(0, 0.3) for xi in x]
    assert paired_permutation_test(x, y, n_permutations=2000)["p_value"] < 0.01


def test_paired_no_effect():
    rng = random.Random(9)
    x = [rng.gauss(0, 1) for _ in range(30)]
    y = [xi + rng.gauss(0, 0.3) for xi in x]
    assert paired_permutation_test(x, y, n_permutations=2000)["p_value"] > 0.05


def test_validation():
    with pytest.raises(ValueError):
        permutation_test([], [1, 2])
    with pytest.raises(ValueError):
        paired_permutation_test([1, 2], [1])
