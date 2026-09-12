"""Entropy pooling: scenario reweighting to honor a mean view."""

import random

import pytest

from quantforge import entropy_pooling_mean, relative_entropy


def _scenarios(n, seed):
    rng = random.Random(seed)
    return [rng.gauss(0, 1) for _ in range(n)]


def test_posterior_sums_to_one_and_hits_target():
    x = _scenarios(500, 1)
    p = entropy_pooling_mean(x, 0.5)
    assert abs(sum(p) - 1.0) < 1e-12
    assert abs(sum(x[i] * p[i] for i in range(len(x))) - 0.5) < 1e-9


def test_no_view_returns_prior():
    x = _scenarios(500, 1)
    prior_mean = sum(x) / len(x)
    p = entropy_pooling_mean(x, prior_mean)
    assert max(abs(pi - 1.0 / len(x)) for pi in p) < 1e-6


def test_relative_entropy_positive_under_view_zero_without():
    x = _scenarios(500, 1)
    q = [1.0 / len(x)] * len(x)
    p_view = entropy_pooling_mean(x, 0.5)
    p_none = entropy_pooling_mean(x, sum(x) / len(x))
    assert relative_entropy(p_view, q) > 0.0
    assert relative_entropy(p_none, q) < 1e-9


def test_stronger_view_costs_more_entropy():
    x = _scenarios(500, 1)
    q = [1.0 / len(x)] * len(x)
    mild = relative_entropy(entropy_pooling_mean(x, 0.3), q)
    strong = relative_entropy(entropy_pooling_mean(x, 0.9), q)
    assert strong > mild


def test_kl_zero_at_equality():
    q = [0.25, 0.25, 0.25, 0.25]
    assert abs(relative_entropy(q, q)) < 1e-15


def test_validation():
    x = _scenarios(100, 1)
    with pytest.raises(ValueError):
        entropy_pooling_mean(x, 10.0)          # target outside scenario range
    with pytest.raises(ValueError):
        entropy_pooling_mean([1.0], 0.5)       # < 2 scenarios
    with pytest.raises(ValueError):
        relative_entropy([0.5, 0.5], [1.0])    # length mismatch
