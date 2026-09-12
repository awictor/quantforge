"""kth-to-default basket probabilities (one-factor Gaussian copula)."""

import pytest

from quantforge import basket_default_distribution, kth_to_default_probability


def test_distribution_sums_to_one():
    d = basket_default_distribution(10, 0.05, 0.2)
    assert abs(sum(d) - 1.0) < 1e-9


def test_monotone_decreasing_in_k():
    p1 = kth_to_default_probability(10, 1, 0.05, 0.2)
    p2 = kth_to_default_probability(10, 2, 0.05, 0.2)
    p5 = kth_to_default_probability(10, 5, 0.05, 0.2)
    assert p1 > p2 > p5


def test_independence_matches_closed_form():
    p1 = kth_to_default_probability(10, 1, 0.05, 0.0)
    assert abs(p1 - (1 - 0.95 ** 10)) < 1e-3


def test_correlation_clusters_defaults():
    # Senior (high-k) trigger rises with correlation; first-to-default falls.
    assert (kth_to_default_probability(20, 5, 0.1, 0.6)
            > kth_to_default_probability(20, 5, 0.1, 0.1))
    assert (kth_to_default_probability(20, 1, 0.1, 0.6)
            < kth_to_default_probability(20, 1, 0.1, 0.1))


def test_expected_defaults_equals_n_times_pd():
    d = basket_default_distribution(50, 0.1, 0.3)
    exp = sum(k * d[k] for k in range(51))
    assert abs(exp - 5.0) < 0.1


def test_validation():
    with pytest.raises(ValueError):
        kth_to_default_probability(10, 0, 0.05, 0.2)      # k < 1
    with pytest.raises(ValueError):
        basket_default_distribution(10, 0.05, 1.0)        # rho = 1
