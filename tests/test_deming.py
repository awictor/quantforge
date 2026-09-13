"""Deming (errors-in-variables) regression."""

import math
import random

import pytest

from quantforge import deming_regression, orthogonal_regression


def test_perfect_line_all_lambdas():
    for lam in (0.5, 1.0, 2.0, 100.0):
        s, b = deming_regression([1, 2, 3, 4, 5], [3, 5, 7, 9, 11], lam)  # y = 2x + 1
        assert abs(s - 2.0) < 1e-9
        assert abs(b - 1.0) < 1e-9


def test_orthogonal_symmetry():
    rng = random.Random(3)
    x = [rng.gauss(0, 1) for _ in range(200)]
    y = [2 * xi + rng.gauss(0, 0.3) for xi in x]
    s_xy, _ = orthogonal_regression(x, y)
    s_yx, _ = orthogonal_regression(y, x)
    assert abs(s_xy * s_yx - 1.0) < 1e-9        # perpendicular fit is symmetric


def test_lambda_infinity_is_ols():
    rng = random.Random(5)
    x = [rng.gauss(0, 1) for _ in range(300)]
    y = [1.5 * xi + 2 + rng.gauss(0, 0.4) for xi in x]
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    ols = sum((x[i] - mx) * (y[i] - my) for i in range(n)) / sum((xi - mx) ** 2 for xi in x)
    big, _ = deming_regression(x, y, lam=1e12)
    assert abs(big - ols) < 1e-4


def test_orthogonal_equals_pca_eigenvector():
    rng = random.Random(5)
    x = [rng.gauss(0, 1) for _ in range(300)]
    y = [1.5 * xi + 2 + rng.gauss(0, 0.4) for xi in x]
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((xi - mx) ** 2 for xi in x) / (n - 1)
    syy = sum((yi - my) ** 2 for yi in y) / (n - 1)
    sxy = sum((x[i] - mx) * (y[i] - my) for i in range(n)) / (n - 1)
    tr, det = sxx + syy, sxx * syy - sxy * sxy
    eig = (tr + math.sqrt(tr * tr - 4 * det)) / 2
    pca_slope = (eig - sxx) / sxy
    o, _ = orthogonal_regression(x, y)
    assert abs(o - pca_slope) < 1e-9


def test_intercept_passes_through_means():
    rng = random.Random(7)
    x = [rng.gauss(0, 1) for _ in range(50)]
    y = [3 * xi - 1 + rng.gauss(0, 0.2) for xi in x]
    s, b = deming_regression(x, y, lam=1.0)
    mx, my = sum(x) / 50, sum(y) / 50
    assert abs(my - (s * mx + b)) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        deming_regression([1, 2], [1])
    with pytest.raises(ValueError):
        deming_regression([1], [1])
    with pytest.raises(ValueError):
        deming_regression([1, 2, 3], [1, 2, 3], lam=0)
    with pytest.raises(ValueError):
        deming_regression([1, 1, 1], [1, 2, 3])       # Sxy = 0
