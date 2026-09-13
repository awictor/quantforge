"""Empirical CDF, sample quantiles, Q-Q pairing."""

import random
import statistics

import pytest

from quantforge import ecdf, quantile, qq_points


def test_ecdf_bounds_and_monotone():
    data = [3, 1, 4, 1, 5, 9, 2, 6]
    assert ecdf(data, 0) == 0.0
    assert ecdf(data, 9) == 1.0
    vals = ecdf(data, [i * 0.5 for i in range(-2, 22)])
    assert all(vals[i] <= vals[i + 1] for i in range(len(vals) - 1))


def test_quantile_linear_matches_statistics():
    data = list(range(1, 101))
    assert quantile(data, 0.5) == 50.5
    q = statistics.quantiles(data, n=4, method="inclusive")
    assert abs(quantile(data, 0.25) - q[0]) < 1e-9
    assert abs(quantile(data, 0.75) - q[2]) < 1e-9


def test_quantile_endpoints():
    data = list(range(1, 101))
    assert quantile(data, 0.0) == 1
    assert quantile(data, 1.0) == 100


def test_quantile_methods():
    d = [10, 20, 30, 40]
    assert quantile(d, 0.5, "lower") == 20
    assert quantile(d, 0.5, "higher") == 30
    assert quantile(d, 0.5, "linear") == 25.0


def test_quantile_recovers_normal():
    rng = random.Random(3)
    big = [rng.gauss(0, 1) for _ in range(10000)]
    assert abs(quantile(big, 0.5)) < 0.1
    assert abs(quantile(big, 0.975) - 1.96) < 0.15


def test_qq_same_distribution_on_diagonal():
    rng = random.Random(5)
    a = [rng.gauss(0, 1) for _ in range(500)]
    b = [rng.gauss(0, 1) for _ in range(500)]
    pts = qq_points(a, b)
    assert max(abs(p[0] - p[1]) for p in pts) < 0.5


def test_qq_scale_difference_slope():
    rng = random.Random(7)
    a = [rng.gauss(0, 1) for _ in range(500)]
    b = [rng.gauss(0, 3) for _ in range(500)]
    pts = qq_points(b, a)                 # b has 3x scale
    far = pts[-50]
    assert abs(far[1] / far[0] - 3) < 0.5


def test_validation():
    with pytest.raises(ValueError):
        ecdf([], 0.0)
    with pytest.raises(ValueError):
        quantile([1, 2], 1.5)
    with pytest.raises(ValueError):
        quantile([1, 2], 0.5, method="bad")
