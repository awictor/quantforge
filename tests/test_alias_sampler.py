"""Walker's alias method for O(1) categorical sampling."""

from collections import Counter

import pytest

from quantforge import AliasSampler


def test_frequencies_match_weights():
    weights = [0.1, 0.2, 0.3, 0.4]
    c = Counter(AliasSampler(weights, seed=42).sample_many(200000))
    for i, w in enumerate(weights):
        assert abs(c[i] / 200000 - w) < 0.01


def test_unnormalized_and_uniform():
    c = Counter(AliasSampler([1, 2, 3, 4], seed=1).sample_many(200000))
    for i, w in enumerate([1, 2, 3, 4]):
        assert abs(c[i] / 200000 - w / 10) < 0.01
    cu = Counter(AliasSampler([1] * 5, seed=7).sample_many(200000))
    assert all(abs(cu[i] / 200000 - 0.2) < 0.01 for i in range(5))


def test_degenerate():
    assert all(x == 2 for x in AliasSampler([0, 0, 5, 0], seed=3).sample_many(1000))


def test_reproducible():
    assert AliasSampler([1, 2, 3], seed=99).sample_many(100) == \
        AliasSampler([1, 2, 3], seed=99).sample_many(100)


def test_validation():
    with pytest.raises(ValueError):
        AliasSampler([])
    with pytest.raises(ValueError):
        AliasSampler([0, 0])
    with pytest.raises(ValueError):
        AliasSampler([-1, 2])
