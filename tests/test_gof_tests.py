"""Goodness-of-fit tests: Jarque-Bera p-value and two-sample KS."""

import random

import pytest

from quantforge import jarque_bera_test, ks_two_sample


def test_normal_sample_not_rejected():
    rng = random.Random(1)
    x = [rng.gauss(0, 1) for _ in range(2000)]
    _, p = jarque_bera_test(x)
    assert p > 0.05


def test_heavy_tail_rejected():
    rng = random.Random(2)
    x = [rng.gauss(0, 1) * (3 if rng.random() < 0.1 else 1) for _ in range(2000)]
    stat, p = jarque_bera_test(x)
    assert stat > 5.99
    assert p < 0.05


def test_ks_same_distribution_not_rejected():
    rng = random.Random(3)
    a = [rng.gauss(0, 1) for _ in range(500)]
    b = [rng.gauss(0, 1) for _ in range(500)]
    d, p = ks_two_sample(a, b)
    assert 0.0 <= d <= 1.0
    assert p > 0.05


def test_ks_location_shift_rejected():
    rng = random.Random(3)
    a = [rng.gauss(0, 1) for _ in range(500)]
    c = [rng.gauss(2, 1) for _ in range(500)]
    d, p = ks_two_sample(a, c)
    assert p < 0.01


def test_ks_scale_difference_rejected():
    rng = random.Random(3)
    a = [rng.gauss(0, 1) for _ in range(500)]
    e = [rng.gauss(0, 3) for _ in range(500)]
    d, p = ks_two_sample(a, e)
    assert p < 0.01


def test_validation():
    with pytest.raises(ValueError):
        ks_two_sample([], [1.0])
    with pytest.raises(ValueError):
        ks_two_sample([1.0], [])
