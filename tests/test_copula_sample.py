"""Gaussian-copula sampler."""

import pytest

from quantforge import gaussian_copula_sample, inverse_transform
from quantforge.mathfns import norm_ppf
from quantforge.copula_stats import spearman_rho


def test_rank_correlation_matches_target():
    s = gaussian_copula_sample([[1.0, 0.7], [0.7, 1.0]], 5000, seed=42)
    x = [r[0] for r in s]
    y = [r[1] for r in s]
    assert abs(spearman_rho(x, y) - 0.7) < 0.05


def test_margins_uniform():
    s = gaussian_copula_sample([[1.0, 0.7], [0.7, 1.0]], 5000, seed=42)
    x = [r[0] for r in s]
    assert abs(sum(x) / len(x) - 0.5) < 0.02
    assert all(0.0 < v < 1.0 for v in x)


def test_independence_near_zero():
    s = gaussian_copula_sample([[1.0, 0.0], [0.0, 1.0]], 5000, seed=7)
    x = [r[0] for r in s]
    y = [r[1] for r in s]
    assert abs(spearman_rho(x, y)) < 0.05


def test_negative_correlation():
    s = gaussian_copula_sample([[1.0, -0.6], [-0.6, 1.0]], 5000, seed=9)
    x = [r[0] for r in s]
    y = [r[1] for r in s]
    assert spearman_rho(x, y) < -0.5


def test_three_variable_shape():
    R = [[1, 0.5, 0.3], [0.5, 1, 0.4], [0.3, 0.4, 1]]
    s = gaussian_copula_sample(R, 3000, seed=1)
    assert len(s) == 3000
    assert all(len(r) == 3 for r in s)


def test_inverse_transform_length():
    s = gaussian_copula_sample([[1.0, 0.5], [0.5, 1.0]], 100, seed=1)
    z = inverse_transform([r[0] for r in s], norm_ppf)
    assert len(z) == 100


def test_validation():
    with pytest.raises(ValueError):
        gaussian_copula_sample([[1.0, 1.5], [1.5, 1.0]], 10)   # not positive-definite
    with pytest.raises(ValueError):
        gaussian_copula_sample([[1.0, 0.5], [0.5, 1.0]], 0)    # n < 1
