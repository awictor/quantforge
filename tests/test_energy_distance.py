"""Energy distance two-sample test."""

import random

import pytest

from quantforge import energy_distance, energy_test


def _ed_ref(a, b):
    n, m = len(a), len(b)
    A = sum(abs(xi - yj) for xi in a for yj in b) / (n * m)
    B = sum(abs(a[i] - a[j]) for i in range(n) for j in range(n)) / (n * n)
    C = sum(abs(b[i] - b[j]) for i in range(m) for j in range(m)) / (m * m)
    return max(2 * A - B - C, 0.0)


def test_identical_is_zero():
    a = [1.0, 2.0, 3.0, 4.0, 5.0]
    assert energy_distance(a, a) == 0.0


def test_symmetric():
    rng = random.Random(1)
    x = [rng.gauss(0, 1) for _ in range(40)]
    y = [rng.gauss(2, 1) for _ in range(30)]
    assert abs(energy_distance(x, y) - energy_distance(y, x)) < 1e-12


def test_matches_reference():
    rng = random.Random(3)
    for _ in range(300):
        n = rng.randint(1, 15)
        m = rng.randint(1, 15)
        a = [rng.gauss(0, 2) for _ in range(n)]
        b = [rng.gauss(0.5, 1) for _ in range(m)]
        assert abs(energy_distance(a, b) - _ed_ref(a, b)) < 1e-9


def test_grows_with_separation():
    rng = random.Random(5)
    a = [rng.gauss(0, 1) for _ in range(200)]
    close = [rng.gauss(0, 1) for _ in range(200)]
    far = [rng.gauss(3, 1) for _ in range(200)]
    assert energy_distance(a, close) < energy_distance(a, far)


def test_permutation_same_vs_different():
    rng = random.Random(7)
    a = [rng.gauss(0, 1) for _ in range(80)]
    b_same = [rng.gauss(0, 1) for _ in range(80)]
    b_diff = [rng.gauss(1.5, 1) for _ in range(80)]
    assert energy_test(a, b_same, n_permutations=299)["p_value"] > 0.05
    assert energy_test(a, b_diff, n_permutations=299)["p_value"] < 0.05


def test_detects_variance_difference():
    # Same mean, different spread: a t-test is blind, energy distance is not.
    rng = random.Random(11)
    a = [rng.gauss(0, 1) for _ in range(150)]
    b = [rng.gauss(0, 3) for _ in range(150)]
    assert energy_test(a, b, n_permutations=299)["p_value"] < 0.05


def test_p_value_bounds():
    rng = random.Random(2)
    a = [rng.gauss(0, 1) for _ in range(30)]
    b = [rng.gauss(0, 1) for _ in range(30)]
    r = energy_test(a, b, n_permutations=199)
    assert 0.0 < r["p_value"] <= 1.0
    assert r["n_permutations"] == 199


def test_validation():
    with pytest.raises(ValueError):
        energy_distance([], [1, 2])
    with pytest.raises(ValueError):
        energy_test([1, 2], [])
