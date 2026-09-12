"""Multivariate Kelly allocation."""

import pytest

from quantforge import (
    kelly_fractions_multivariate, kelly_growth_rate_multivariate,
    kelly_fraction_continuous,
)


def test_scalar_reduction():
    f = kelly_fractions_multivariate([0.08], [[0.04]])
    assert abs(f[0] - kelly_fraction_continuous(0.08, 0.04)) < 1e-12


def test_diagonal_equals_per_asset_kelly():
    mu = [0.08, 0.05, 0.10]
    var = [0.04, 0.02, 0.09]
    cov = [[var[i] if i == j else 0.0 for j in range(3)] for i in range(3)]
    f = kelly_fractions_multivariate(mu, cov)
    for i in range(3):
        assert abs(f[i] - mu[i] / var[i]) < 1e-12


def test_growth_maximized_at_optimum():
    mu = [0.08, 0.05]
    cov = [[0.04, 0.012], [0.012, 0.02]]
    fstar = kelly_fractions_multivariate(mu, cov)
    g0 = kelly_growth_rate_multivariate(mu, cov, fstar)
    for d0 in (-0.5, 0.0, 0.5):
        for d1 in (-0.5, 0.0, 0.5):
            g = kelly_growth_rate_multivariate(mu, cov, [fstar[0] + d0, fstar[1] + d1])
            assert g <= g0 + 1e-12


def test_fractional_scaling():
    mu = [0.08, 0.05]
    cov = [[0.04, 0.012], [0.012, 0.02]]
    full = kelly_fractions_multivariate(mu, cov)
    half = kelly_fractions_multivariate(mu, cov, fraction=0.5)
    for i in range(2):
        assert abs(half[i] - 0.5 * full[i]) < 1e-12


def test_positive_correlation_reduces_total_leverage():
    mu = [0.06, 0.06]
    corr = [[0.04, 0.036], [0.036, 0.04]]   # correlation 0.9
    indep = [[0.04, 0.0], [0.0, 0.04]]
    assert sum(kelly_fractions_multivariate(mu, corr)) < sum(
        kelly_fractions_multivariate(mu, indep))


def test_validation():
    with pytest.raises(ValueError):
        kelly_fractions_multivariate([0.1, 0.1], [[0.04]])       # shape mismatch
    with pytest.raises(ValueError):
        kelly_fractions_multivariate([0.1, 0.1], [[0.0, 0.0], [0.0, 0.0]])  # singular
    with pytest.raises(ValueError):
        kelly_growth_rate_multivariate([0.1], [[0.04]], [1.0, 2.0])  # leverage mismatch
