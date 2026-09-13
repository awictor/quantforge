"""Compensated summation and stable running statistics."""

import math
import random
import statistics

import pytest

from quantforge import kahan_sum, neumaier_sum, accurate_dot, welford


def test_neumaier_handles_catastrophic_cancellation():
    assert neumaier_sum([1.0, 1e100, 1.0, -1e100]) == 2.0


def test_matches_fsum():
    rng = random.Random(3)
    for _ in range(50):
        vals = [rng.uniform(-1, 1) * 10 ** rng.randint(-8, 8) for _ in range(1000)]
        target = math.fsum(vals)
        assert abs(neumaier_sum(vals) - target) < 1e-6 * max(1.0, abs(target))


def test_kahan_beats_naive():
    n = 1000000
    naive = 0.0
    for _ in range(n):
        naive += 0.1
    k = kahan_sum([0.1] * n)
    assert abs(k - 100000) < abs(naive - 100000)
    assert abs(k - 100000) < 1e-6


def test_accurate_dot():
    assert accurate_dot([1e8, 1, -1e8], [1, 1, 1]) == 1.0


def test_welford_matches_statistics():
    rng = random.Random(5)
    data = [rng.gauss(5, 2) for _ in range(1000)]
    m, v, n = welford(data)
    assert abs(m - statistics.mean(data)) < 1e-9
    assert abs(v - statistics.variance(data)) < 1e-9
    assert n == 1000


def test_welford_stable_for_large_mean():
    rng = random.Random(7)
    big = [1e9 + rng.gauss(0, 1) for _ in range(1000)]
    _, v, _ = welford(big)
    # Naive E[x^2]-E[x]^2 collapses to ~0 here; Welford stays near the true variance 1.
    assert 0.5 < v < 2.0


def test_validation():
    with pytest.raises(ValueError):
        accurate_dot([1, 2], [1])
    with pytest.raises(ValueError):
        welford([])
