"""Tests for weighted reservoir and with-replacement sampling."""

import math
from collections import Counter

import pytest

from quantforge.weighted_reservoir import (
    weighted_reservoir_sample,
    weighted_sample_with_replacement,
)


def test_k_one_inclusion_proportional_to_weight():
    # Efraimidis-Spirakis with k=1 reduces to weighted choice: P(pick i) = w_i / sum(w).
    items = ["a", "b", "c", "d"]
    weights = [1, 2, 3, 4]
    counts = Counter()
    trials = 30000
    for s in range(trials):
        (pick,) = weighted_reservoir_sample(items, weights, 1, seed=s)
        counts[pick] += 1
    total = sum(weights)
    for it, w in zip(items, weights):
        assert abs(counts[it] / trials - w / total) < 0.02


def test_returns_k_distinct_items():
    items = list(range(10))
    weights = [i + 1 for i in range(10)]
    for s in range(500):
        r = weighted_reservoir_sample(items, weights, 4, seed=s)
        assert len(r) == 4
        assert len(set(r)) == 4
        assert all(x in items for x in r)


def test_k_geq_positive_count_returns_all():
    items = ["a", "b", "c"]
    weights = [5, 5, 5]
    r = weighted_reservoir_sample(items, weights, 10, seed=1)
    assert set(r) == set(items)


def test_zero_weight_items_never_selected():
    items = ["x", "y", "z"]
    weights = [0.0, 5.0, 0.0]
    for s in range(200):
        r = weighted_reservoir_sample(items, weights, 2, seed=s)
        assert r == ["y"]  # only one positive-weight item exists


def test_k_zero_returns_empty():
    assert weighted_reservoir_sample(["a", "b"], [1, 1], 0) == []


def test_reservoir_deterministic_for_seed():
    items = list(range(20))
    weights = [1.0] * 20
    a = weighted_reservoir_sample(items, weights, 5, seed=99)
    b = weighted_reservoir_sample(items, weights, 5, seed=99)
    assert a == b


def test_reservoir_length_mismatch_raises():
    with pytest.raises(ValueError):
        weighted_reservoir_sample(["a"], [1, 2], 1)


def test_reservoir_negative_k_raises():
    with pytest.raises(ValueError):
        weighted_reservoir_sample(["a", "b"], [1, 1], -1)


def test_with_replacement_frequency_matches_weights():
    items = ["a", "b", "c", "d"]
    weights = [1, 2, 3, 4]
    n = 200000
    counts = Counter(weighted_sample_with_replacement(items, weights, n, seed=1))
    total = sum(weights)
    for it, w in zip(items, weights):
        assert abs(counts[it] / n - w / total) < 0.01


def test_with_replacement_length_and_membership():
    items = ["a", "b"]
    r = weighted_sample_with_replacement(items, [1, 1], 50, seed=7)
    assert len(r) == 50
    assert all(x in items for x in r)


def test_with_replacement_can_repeat():
    # single positive weight -> every draw is that item
    r = weighted_sample_with_replacement(["a", "b"], [1.0, 0.0], 10, seed=3)
    assert r == ["a"] * 10


def test_with_replacement_deterministic_for_seed():
    items = list(range(5))
    weights = [1, 2, 3, 4, 5]
    a = weighted_sample_with_replacement(items, weights, 30, seed=42)
    b = weighted_sample_with_replacement(items, weights, 30, seed=42)
    assert a == b


def test_with_replacement_negative_weight_raises():
    with pytest.raises(ValueError):
        weighted_sample_with_replacement(["a", "b"], [1, -1], 5)


def test_with_replacement_zero_total_raises():
    with pytest.raises(ValueError):
        weighted_sample_with_replacement(["a", "b"], [0, 0], 5)


def test_with_replacement_length_mismatch_raises():
    with pytest.raises(ValueError):
        weighted_sample_with_replacement(["a"], [1, 2], 3)


def test_two_item_reservoir_inclusion_probability():
    # For k=1, two items: closed-form P(a) = wa/(wa+wb). Check wa=3, wb=1 -> 0.75.
    counts = Counter()
    trials = 20000
    for s in range(trials):
        (pick,) = weighted_reservoir_sample(["a", "b"], [3.0, 1.0], 1, seed=s)
        counts[pick] += 1
    assert abs(counts["a"] / trials - 0.75) < 0.02
