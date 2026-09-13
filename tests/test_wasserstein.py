"""1-D Wasserstein (earth-mover) distance."""

import random

import pytest

from quantforge import wasserstein_distance, wasserstein1_sorted


def test_constant_shift():
    x = [1, 2, 3, 4, 5]
    y = [xi + 3 for xi in x]
    assert abs(wasserstein1_sorted(x, y) - 3.0) < 1e-12
    assert abs(wasserstein_distance(x, y) - 3.0) < 1e-9


def test_identical_zero():
    x = [1, 2, 3, 4, 5]
    assert wasserstein_distance(x, x) == 0.0


def test_sorted_matches_general_equal_length():
    rng = random.Random(3)
    a = [rng.gauss(0, 1) for _ in range(100)]
    b = [rng.gauss(1, 1) for _ in range(100)]
    assert abs(wasserstein1_sorted(a, b) - wasserstein_distance(a, b)) < 1e-9


def test_point_masses():
    assert wasserstein_distance([0.0], [1.0]) == 1.0


def test_symmetric():
    rng = random.Random(5)
    a = [rng.gauss(0, 1) for _ in range(80)]
    b = [rng.gauss(2, 1) for _ in range(60)]
    assert abs(wasserstein_distance(a, b) - wasserstein_distance(b, a)) < 1e-12


def test_w2_geq_w1():
    rng = random.Random(7)
    a = [rng.gauss(0, 1) for _ in range(100)]
    b = [rng.gauss(1, 2) for _ in range(100)]
    assert wasserstein_distance(a, b, p=2) >= wasserstein_distance(a, b, p=1) - 1e-12


def test_uniform_shift_large_sample():
    rng = random.Random(9)
    a = [rng.random() for _ in range(5000)]
    b = [rng.random() + 2 for _ in range(5000)]
    assert abs(wasserstein_distance(a, b) - 2.0) < 0.05


def test_validation():
    with pytest.raises(ValueError):
        wasserstein1_sorted([1, 2], [1])
    with pytest.raises(ValueError):
        wasserstein_distance([], [1])
    with pytest.raises(ValueError):
        wasserstein_distance([1], [1], p=0)
