"""Financial turbulence and absorption ratio."""

import random

import pytest

from quantforge import turbulence, turbulence_series, absorption_ratio


def test_turbulence_mean_equals_dimensions():
    rng = random.Random(1)
    n = 4
    data = [[rng.gauss(0, 1) for _ in range(n)] for _ in range(3000)]
    ts = turbulence_series(data)
    assert abs(sum(ts) / len(ts) - n) < 0.3


def test_outlier_high_turbulence():
    mean = [0.0] * 4
    cov = [[1.0 if i == j else 0.0 for j in range(4)] for i in range(4)]
    normal = turbulence([0.1, 0.1, 0.1, 0.1], mean, cov)
    outlier = turbulence([5, -5, 5, -5], mean, cov)
    assert outlier > 50 * normal


def test_turbulence_non_negative():
    rng = random.Random(2)
    data = [[rng.gauss(0, 1) for _ in range(3)] for _ in range(200)]
    assert all(t >= 0.0 for t in turbulence_series(data))


def test_absorption_ratio_even():
    cov = [[1.0 if i == j else 0.0 for j in range(4)] for i in range(4)]
    assert abs(absorption_ratio(cov, 1) - 0.25) < 1e-9


def test_absorption_ratio_dominant_near_one():
    cov = [[100.0 if i == j == 0 else (1.0 if i == j else 0.0) for j in range(4)]
           for i in range(4)]
    assert absorption_ratio(cov, 1) > 0.9


def test_absorption_ratio_in_unit_interval():
    rng = random.Random(3)
    A = [[rng.gauss(0, 1) for _ in range(4)] for _ in range(4)]
    cov = [[sum(A[i][k] * A[j][k] for k in range(4)) for j in range(4)]
           for i in range(4)]
    ar = absorption_ratio(cov, 2)
    assert 0.0 <= ar <= 1.0


def test_validation():
    with pytest.raises(ValueError):
        turbulence([1.0, 2.0], [0.0], [[1.0]])
    with pytest.raises(ValueError):
        absorption_ratio([[1.0, 0.0], [0.0, 1.0]], n_factors=3)
    with pytest.raises(ValueError):
        turbulence_series([[1.0, 2.0]])
