"""LOWESS locally-weighted scatterplot smoothing."""

import math
import random

import pytest

from quantforge import lowess


def test_line_reproduced():
    x = [float(i) for i in range(50)]
    y = [2 * xi + 1 for xi in x]
    sm = lowess(x, y, frac=0.5)
    assert max(abs(sm[i] - y[i]) for i in range(50)) < 1e-6


def test_smooths_noisy_curve():
    rng = random.Random(3)
    x = [i * 0.1 for i in range(100)]
    true = [math.sin(xi) for xi in x]
    yn = [true[i] + rng.gauss(0, 0.2) for i in range(100)]
    sm = lowess(x, yn, frac=0.3)
    raw = sum((yn[i] - true[i]) ** 2 for i in range(100)) / 100
    smoothed = sum((sm[i] - true[i]) ** 2 for i in range(100)) / 100
    assert smoothed < raw


def test_robust_to_outlier():
    x = [float(i) for i in range(50)]
    y = [0.5 * xi for xi in x]
    y[25] += 100
    sm = lowess(x, y, frac=0.4, iterations=3)
    assert abs(sm[24] - 12.0) < 0.5
    assert abs(sm[26] - 13.0) < 0.5


def test_larger_frac_is_smoother():
    rng = random.Random(5)
    x = [float(i) for i in range(50)]
    yn = [math.sin(i * 0.2) + rng.gauss(0, 0.2) for i in range(50)]
    sm_s = lowess(x, yn, frac=0.1)
    sm_b = lowess(x, yn, frac=0.8)
    tv = lambda v: sum(abs(v[i + 1] - v[i]) for i in range(len(v) - 1))
    assert tv(sm_b) < tv(sm_s)


def test_preserves_input_order():
    sm = lowess([3.0, 1.0, 2.0], [6.0, 2.0, 4.0], frac=1.0)
    assert abs(sm[0] - 6.0) < 1e-6 and abs(sm[1] - 2.0) < 1e-6 and abs(sm[2] - 4.0) < 1e-6


def test_validation():
    with pytest.raises(ValueError):
        lowess([1.0], [1.0])
    with pytest.raises(ValueError):
        lowess([1.0, 2.0], [1.0, 2.0], frac=1.5)
