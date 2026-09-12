"""Rank-based dependence: Kendall's tau, Spearman's rho, empirical copula."""

import math
import random

import pytest

from quantforge import kendall_tau, spearman_rho, pseudo_observations


def test_perfect_increasing():
    x = list(range(20))
    y = [2 * v + 3 for v in x]
    assert abs(kendall_tau(x, y) - 1.0) < 1e-12
    assert abs(spearman_rho(x, y) - 1.0) < 1e-12


def test_perfect_decreasing():
    x = list(range(20))
    y = [-v for v in x]
    assert abs(kendall_tau(x, y) + 1.0) < 1e-12
    assert abs(spearman_rho(x, y) + 1.0) < 1e-12


def test_spearman_monotone_transform_invariant():
    x = list(range(20))
    y = [math.exp(v / 5.0) for v in x]     # strictly increasing transform
    assert abs(spearman_rho(x, y) - 1.0) < 1e-12


def test_independence_near_zero():
    rng = random.Random(1)
    a = [rng.gauss(0, 1) for _ in range(2000)]
    b = [rng.gauss(0, 1) for _ in range(2000)]
    assert abs(kendall_tau(a, b)) < 0.06
    assert abs(spearman_rho(a, b)) < 0.06


def test_pseudo_observations_in_open_unit_interval():
    p = pseudo_observations([50, 10, 30, 20, 40])
    assert all(0.0 < v < 1.0 for v in p)
    assert pseudo_observations([1, 2, 3, 4]) == [0.2, 0.4, 0.6, 0.8]


def test_pseudo_observations_average_ties():
    # [1,1,2] -> ranks 1.5,1.5,3 -> /(n+1)=/4
    assert pseudo_observations([1, 1, 2]) == [0.375, 0.375, 0.75]


def test_validation():
    with pytest.raises(ValueError):
        kendall_tau([1.0], [1.0])           # < 2 points
    with pytest.raises(ValueError):
        spearman_rho([1.0, 2.0], [1.0])     # length mismatch
    with pytest.raises(ValueError):
        spearman_rho([5.0] * 4, [1, 2, 3, 4])  # no rank variation in x
