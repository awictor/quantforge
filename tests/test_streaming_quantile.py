"""Streaming quantile (P-square) and reservoir sampling."""

import random

import pytest

from quantforge import P2Quantile, reservoir_sample


def test_p2_median_uniform():
    rng = random.Random(3)
    q = P2Quantile(0.5)
    for _ in range(100000):
        q.update(rng.random())
    assert abs(q.value() - 0.5) < 0.02


def test_p2_matches_true_percentile():
    rng = random.Random(3)
    data = [rng.random() for _ in range(100000)]
    q = P2Quantile(0.95)
    for x in data:
        q.update(x)
    true95 = sorted(data)[int(0.95 * len(data))]
    assert abs(q.value() - true95) < 0.01


def test_p2_normal_p90():
    rng = random.Random(5)
    q = P2Quantile(0.9)
    for _ in range(100000):
        q.update(rng.gauss(0, 1))
    assert abs(q.value() - 1.2816) < 0.05


def test_p2_exact_below_five():
    q = P2Quantile(0.5)
    for x in [3, 1, 2]:
        q.update(x)
    assert q.value() == 2


def test_reservoir_uniform_coverage():
    counts = [0] * 100
    for seed in range(2000):
        for x in reservoir_sample(range(100), 10, seed=seed):
            counts[x] += 1
    # Expected ~200 per item; uniform sampler keeps the spread tight.
    assert min(counts) > 120
    assert max(counts) < 280


def test_reservoir_small_stream_returns_all():
    assert sorted(reservoir_sample([1, 2, 3], 10)) == [1, 2, 3]


def test_reservoir_deterministic():
    assert reservoir_sample(range(1000), 5, seed=42) == reservoir_sample(range(1000), 5, seed=42)


def test_validation():
    with pytest.raises(ValueError):
        P2Quantile(0.0)
    with pytest.raises(ValueError):
        reservoir_sample([1, 2, 3], 0)
