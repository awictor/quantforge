"""RANSAC robust line fitting."""

import random

import pytest

from quantforge import ransac_line


def test_recovers_line_with_60pct_outliers():
    rng = random.Random(3)
    x = [float(i) for i in range(100)]
    y = [2 * x[i] + 5 + rng.gauss(0, 0.5) for i in range(100)]
    for i in random.Random(1).sample(range(100), 60):
        y[i] = rng.uniform(-200, 400)
    r = ransac_line(x, y, threshold=3.0, n_iterations=500)
    assert abs(r["slope"] - 2) < 0.1
    assert abs(r["intercept"] - 5) < 1.0


def test_clean_data_exact():
    x = [float(i) for i in range(20)]
    y = [3 * x[i] - 2 for i in range(20)]
    r = ransac_line(x, y, threshold=0.01)
    assert abs(r["slope"] - 3) < 1e-6 and abs(r["intercept"] + 2) < 1e-6
    assert r["n_inliers"] == 20


def test_deterministic():
    x = [float(i) for i in range(20)]
    y = [3 * x[i] - 2 for i in range(20)]
    a = ransac_line(x, y, threshold=1.0, seed=42)
    b = ransac_line(x, y, threshold=1.0, seed=42)
    assert a["slope"] == b["slope"] and a["intercept"] == b["intercept"]


def test_majority_structure():
    rng = random.Random(7)
    x = [float(i) for i in range(60)]
    y = [x[i] + 10 + rng.gauss(0, 0.2) if i < 40 else rng.uniform(0, 100)
         for i in range(60)]
    r = ransac_line(x, y, threshold=1.0, n_iterations=500)
    assert abs(r["slope"] - 1) < 0.1


def test_validation():
    with pytest.raises(ValueError):
        ransac_line([1.0], [1.0], threshold=1.0)
    with pytest.raises(ValueError):
        ransac_line([1.0, 2.0], [1.0, 2.0], threshold=0)
