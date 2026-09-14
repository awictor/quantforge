"""Random variate samplers: normal, exponential, gamma, Poisson."""

import statistics

import pytest

from quantforge import sample_normal, sample_exponential, sample_gamma, sample_poisson

N = 100000


def test_normal_moments():
    x = sample_normal(N, mu=5, sigma=2, seed=1)
    assert abs(statistics.mean(x) - 5) < 0.05
    assert abs(statistics.pstdev(x) - 2) < 0.05


def test_exponential_moments():
    x = sample_exponential(N, rate=0.5, seed=2)
    assert abs(statistics.mean(x) - 2) < 0.05
    assert abs(statistics.pvariance(x) - 4) < 0.3


def test_gamma_moments():
    for shape in (0.5, 2.0, 7.5):
        x = sample_gamma(N, shape, scale=2.0, seed=3)
        assert abs(statistics.mean(x) - shape * 2) < 0.1
        assert abs(statistics.pvariance(x) - shape * 4) < 0.6


def test_poisson_moments():
    x = sample_poisson(N, lam=4.0, seed=4)
    assert abs(statistics.mean(x) - 4) < 0.05
    assert abs(statistics.pvariance(x) - 4) < 0.1


def test_reproducible():
    assert sample_normal(100, seed=9) == sample_normal(100, seed=9)
    assert sample_gamma(100, 2.0, seed=9) == sample_gamma(100, 2.0, seed=9)


def test_validation():
    with pytest.raises(ValueError):
        sample_normal(10, 0, -1)
    with pytest.raises(ValueError):
        sample_exponential(10, 0)
    with pytest.raises(ValueError):
        sample_gamma(10, -1)
    with pytest.raises(ValueError):
        sample_poisson(10, -1)
