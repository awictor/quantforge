"""Empirical tail dependence and exceedance correlation."""

import random

import pytest

from quantforge import (
    upper_tail_dependence, lower_tail_dependence, exceedance_correlation,
)


def test_comonotone_tail_dependence_near_one():
    rng = random.Random(1)
    x = [rng.gauss(0, 1) for _ in range(3000)]
    y = [3 * v + 1 for v in x]                 # strictly increasing
    assert upper_tail_dependence(x, y, 0.9) > 0.98
    assert lower_tail_dependence(x, y, 0.1) > 0.98


def test_independence_matches_null_fraction():
    rng = random.Random(1)
    a = [rng.gauss(0, 1) for _ in range(5000)]
    b = [rng.gauss(0, 1) for _ in range(5000)]
    # Under independence P(U>q | V>q) = P(U>q) = 1 - q.
    assert abs(upper_tail_dependence(a, b, 0.9) - 0.1) < 0.05


def test_more_extreme_threshold_shrinks_independent_tail():
    rng = random.Random(1)
    a = [rng.gauss(0, 1) for _ in range(5000)]
    b = [rng.gauss(0, 1) for _ in range(5000)]
    assert upper_tail_dependence(a, b, 0.99) < upper_tail_dependence(a, b, 0.9) + 0.05


def test_common_shock_has_positive_tail_dependence():
    rng = random.Random(2)
    xs, ys = [], []
    for _ in range(5000):
        shock = rng.gauss(0, 1)
        xs.append(shock + 0.5 * rng.gauss(0, 1))
        ys.append(shock + 0.5 * rng.gauss(0, 1))
    assert upper_tail_dependence(xs, ys, 0.9) > 0.2


def test_exceedance_correlation_high_when_comonotone():
    rng = random.Random(3)
    x = [rng.gauss(0, 1) for _ in range(3000)]
    y = [3 * v + 1 for v in x]
    assert exceedance_correlation(x, y, 0.9, "upper") > 0.99


def test_validation():
    x = [1.0, 2.0, 3.0, 4.0]
    with pytest.raises(ValueError):
        upper_tail_dependence(x, x, 1.5)          # q outside (0,1)
    with pytest.raises(ValueError):
        upper_tail_dependence(x, [1.0], 0.9)      # length mismatch
    with pytest.raises(ValueError):
        exceedance_correlation(x, x, 0.9, "sideways")  # bad tail
