"""Theil-Sen robust regression."""

import random

import pytest

from quantforge import theil_sen


def _ols(x, y):
    n = len(x)
    xm = sum(x) / n
    ym = sum(y) / n
    s = sum((x[i] - xm) * (y[i] - ym) for i in range(n)) / sum((xi - xm) ** 2 for xi in x)
    return s, ym - s * xm


def test_exact_on_noiseless_line():
    x = list(range(20))
    y = [3.0 * xi + 7.0 for xi in x]
    s, b = theil_sen(x, y)
    assert abs(s - 3.0) < 1e-12
    assert abs(b - 7.0) < 1e-12


def test_close_to_ols_on_clean_data():
    rng = random.Random(1)
    x = [rng.gauss(0, 1) for _ in range(200)]
    y = [2.0 * xi + 1.0 + rng.gauss(0, 0.3) for xi in x]
    s, _ = theil_sen(x, y)
    ols_s, _ = _ols(x, y)
    assert abs(s - ols_s) < 0.1
    assert abs(s - 2.0) < 0.1


def test_robust_where_ols_breaks():
    rng = random.Random(2)
    x = [float(i) for i in range(100)]
    y = [1.5 * xi + 2.0 for xi in x]
    yc = list(y)
    for _ in range(20):
        yc[rng.randint(0, 99)] += 500.0 * rng.choice([-1, 1])
    s, _ = theil_sen(x, yc)
    ols_s, _ = _ols(x, yc)
    assert abs(s - 1.5) < 0.2       # Theil-Sen holds the true slope
    assert abs(ols_s - 1.5) > 0.3   # OLS is dragged away


def test_duplicate_x_pairs_skipped():
    s, _ = theil_sen([1, 1, 2, 3], [5, 9, 4, 6])
    assert isinstance(s, float)


def test_validation():
    with pytest.raises(ValueError):
        theil_sen([1.0, 2.0], [1.0])          # length mismatch
    with pytest.raises(ValueError):
        theil_sen([1.0], [2.0])               # < 2 points
    with pytest.raises(ValueError):
        theil_sen([1.0, 1.0, 1.0], [2, 3, 4]) # all x identical
