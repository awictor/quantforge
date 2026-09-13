"""Time-series entropy: approximate, sample, permutation."""

import math
import random

import pytest

from quantforge import approximate_entropy, sample_entropy, permutation_entropy


def _sine(n, period=20):
    return [math.sin(2 * math.pi * i / period) for i in range(n)]


def _random(n, seed):
    rng = random.Random(seed)
    return [rng.gauss(0, 1) for _ in range(n)]


def test_apen_regular_below_random():
    assert approximate_entropy(_sine(400)) < approximate_entropy(_random(400, 1))


def test_sampen_regular_below_random():
    assert sample_entropy(_sine(400)) < sample_entropy(_random(400, 1))


def test_sampen_positive_for_noise():
    assert sample_entropy(_random(400, 2)) > 0.5


def test_permutation_entropy_monotone_is_zero():
    assert permutation_entropy(list(range(200))) == 0.0
    assert permutation_entropy(list(range(200, 0, -1))) == 0.0


def test_permutation_entropy_random_near_one():
    assert permutation_entropy(_random(2000, 3)) > 0.95


def test_permutation_entropy_in_unit_interval():
    for s in range(5):
        pe = permutation_entropy(_random(500, s), m=4)
        assert 0.0 <= pe <= 1.0


def test_permutation_entropy_monotone_transform_invariant():
    x = _random(500, 7)
    pe1 = permutation_entropy(x, m=4)
    pe2 = permutation_entropy([math.exp(v) for v in x], m=4)   # exp is increasing
    assert abs(pe1 - pe2) < 1e-12


def test_validation():
    with pytest.raises(ValueError):
        approximate_entropy([1.0, 2.0])                 # too short
    with pytest.raises(ValueError):
        sample_entropy([1.0, 1.0, 1.0, 1.0], r=0.0)     # zero tolerance
    with pytest.raises(ValueError):
        permutation_entropy([1.0, 2.0], m=3)            # too short
    with pytest.raises(ValueError):
        permutation_entropy(_random(50, 1), m=1)        # m < 2
