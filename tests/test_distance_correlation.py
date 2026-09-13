"""Distance correlation (Szekely-Rizzo)."""

import math
import random

import pytest

from quantforge import distance_correlation, distance_covariance, distance_variance


def _brute_dcor(x, y):
    n = len(x)

    def dc(s):
        D = [[abs(s[i] - s[j]) for j in range(n)] for i in range(n)]
        rm = [sum(D[i]) / n for i in range(n)]
        cm = [sum(D[i][j] for i in range(n)) / n for j in range(n)]
        g = sum(rm) / n
        return [[D[i][j] - rm[i] - cm[j] + g for j in range(n)] for i in range(n)]

    A, B = dc(x), dc(y)
    dcov = sum(A[i][j] * B[i][j] for i in range(n) for j in range(n)) / (n * n)
    dvx = sum(A[i][j] ** 2 for i in range(n) for j in range(n)) / (n * n)
    dvy = sum(B[i][j] ** 2 for i in range(n) for j in range(n)) / (n * n)
    return math.sqrt(max(dcov, 0) / math.sqrt(dvx * dvy)) if dvx * dvy > 0 else 0.0


def test_linear_is_one():
    assert abs(distance_correlation([0, 1, 2, 3, 4], [0, 1, 2, 3, 4]) - 1.0) < 1e-12
    # Scale-invariant: dCor of y = 2x is also 1.
    assert abs(distance_correlation([0, 1, 2, 3, 4], [0, 2, 4, 6, 8]) - 1.0) < 1e-12


def test_independence_is_small():
    rng = random.Random(3)
    x = [rng.gauss(0, 1) for _ in range(400)]
    y = [rng.gauss(0, 1) for _ in range(400)]
    assert distance_correlation(x, y) < 0.15


def test_detects_nonlinear_where_pearson_fails():
    rng = random.Random(9)
    x = [rng.uniform(-3, 3) for _ in range(400)]
    y = [xi * xi for xi in x]
    mx, my = sum(x) / len(x), sum(y) / len(y)
    cov = sum((x[i] - mx) * (y[i] - my) for i in range(len(x)))
    sx = math.sqrt(sum((v - mx) ** 2 for v in x))
    sy = math.sqrt(sum((v - my) ** 2 for v in y))
    pearson = cov / (sx * sy)
    assert abs(pearson) < 0.2                    # Pearson is blind to the parabola
    assert distance_correlation(x, y) > 0.4      # distance correlation sees it


def test_matches_brute_force():
    rng = random.Random(11)
    for _ in range(200):
        n = rng.randint(2, 12)
        x = [rng.gauss(0, 1) for _ in range(n)]
        y = [rng.gauss(0, 1) for _ in range(n)]
        assert abs(distance_correlation(x, y) - _brute_dcor(x, y)) < 1e-9


def test_in_unit_interval():
    rng = random.Random(1)
    x = [rng.gauss(0, 1) for _ in range(50)]
    y = [rng.gauss(0, 1) for _ in range(50)]
    dc = distance_correlation(x, y)
    assert 0.0 <= dc <= 1.0


def test_constant_is_zero():
    assert distance_correlation([5, 5, 5, 5], [1, 2, 3, 4]) == 0.0
    assert distance_variance([5, 5, 5, 5]) == 0.0


def test_dvar_is_dcov_with_self():
    x = [1.0, 3.0, 2.0, 7.0, 4.0]
    assert abs(distance_variance(x) - distance_covariance(x, x)) < 1e-12


def test_validation():
    with pytest.raises(ValueError):
        distance_correlation([1, 2], [1])
    with pytest.raises(ValueError):
        distance_correlation([1], [1])
