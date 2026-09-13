"""Siegel repeated-median regression."""

import random

import pytest

from quantforge import repeated_median_regression, theil_sen


def test_clean_line_exact():
    x = list(range(10))
    y = [2 * xi + 1 for xi in x]
    s, b = repeated_median_regression(x, y)
    assert abs(s - 2.0) < 1e-9 and abs(b - 1.0) < 1e-9


def test_survives_40_percent_outliers():
    x = list(range(20))
    y = [3 * xi - 5 for xi in x]
    for i in [1, 4, 7, 10, 13, 16, 2, 5]:            # 8/20 = 40% corrupted
        y[i] = random.Random(i).uniform(-500, 500)
    s, b = repeated_median_regression(x, y)
    assert abs(s - 3.0) < 1e-6 and abs(b + 5.0) < 1e-6


def test_decreasing_line():
    s, b = repeated_median_regression([1, 2, 3, 4, 5], [10, 8, 6, 4, 2])
    assert abs(s + 2.0) < 1e-9 and abs(b - 12.0) < 1e-9


def test_matches_theil_sen_on_clean():
    rng = random.Random(3)
    x = [rng.uniform(0, 10) for _ in range(30)]
    y = [1.5 * xi + 2 + rng.gauss(0, 0.3) for xi in x]
    s, _ = repeated_median_regression(x, y)
    ts, _ = theil_sen(x, y)
    assert abs(s - ts) < 0.1


def test_single_outlier_unmoved():
    x = list(range(10))
    y = [2 * xi + 1 for xi in x]
    y[5] = 1000
    s, b = repeated_median_regression(x, y)
    assert abs(s - 2.0) < 1e-9 and abs(b - 1.0) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        repeated_median_regression([1], [1])
    with pytest.raises(ValueError):
        repeated_median_regression([5, 5, 5], [1, 2, 3])     # all x identical
