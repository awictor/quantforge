"""Passing-Bablok method-comparison regression."""

import random

import pytest

from quantforge import passing_bablok_regression


def test_perfect_line():
    s, b = passing_bablok_regression([1, 2, 3, 4, 5], [3, 5, 7, 9, 11])
    assert abs(s - 2.0) < 1e-9
    assert abs(b - 1.0) < 1e-9


def test_robust_to_outlier():
    x = list(range(1, 11))
    y = [2 * xi + 1 for xi in x]
    y[5] = 100          # gross outlier
    s, b = passing_bablok_regression(x, y)
    assert abs(s - 2.0) < 1e-9         # median slope unmoved
    assert abs(b - 1.0) < 1e-9


def test_decreasing_line():
    s, b = passing_bablok_regression([1, 2, 3, 4, 5], [10, 8, 6, 4, 2])
    assert abs(s - (-2.0)) < 1e-9
    assert abs(b - 12.0) < 1e-9


def test_recovers_slope_on_noisy_data():
    rng = random.Random(3)
    x = [rng.uniform(1, 10) for _ in range(60)]
    y = [1.5 * xi + 2 + rng.gauss(0, 0.5) for xi in x]
    s, _ = passing_bablok_regression(x, y)
    assert abs(s - 1.5) < 0.15


def test_symmetry():
    rng = random.Random(5)
    x = [rng.uniform(1, 10) for _ in range(50)]
    y = [1.5 * xi + 2 + rng.gauss(0, 0.4) for xi in x]
    s1, _ = passing_bablok_regression(x, y)
    s2, _ = passing_bablok_regression(y, x)
    assert abs(s1 * s2 - 1.0) < 0.1        # reciprocal-slope symmetry


def test_no_index_errors_on_random():
    rng = random.Random(7)
    for _ in range(1000):
        n = rng.randint(2, 15)
        xs = [rng.uniform(-5, 5) for _ in range(n)]
        ys = [rng.uniform(-5, 5) for _ in range(n)]
        try:
            passing_bablok_regression(xs, ys)
        except ValueError:
            pass


def test_validation():
    with pytest.raises(ValueError):
        passing_bablok_regression([1, 2], [1])
    with pytest.raises(ValueError):
        passing_bablok_regression([1], [1])
