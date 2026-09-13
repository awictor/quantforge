"""Generic information criteria (AIC / AICc / BIC / HQIC)."""

import math

import pytest

from quantforge import gaussian_log_likelihood, aic, aicc, bic, hqic


def test_aic_formula():
    assert aic(-100.0, 3) == 206.0


def test_bic_exceeds_aic_for_large_n():
    ll, k = -100.0, 3
    assert bic(ll, k, 8) > aic(ll, k)
    assert bic(ll, k, 100) > aic(ll, k)


def test_aicc_approaches_aic():
    ll, k = -100.0, 3
    near = aicc(ll, k, 20)
    far = aicc(ll, k, 100000)
    assert near > aic(ll, k)
    assert abs(far - aic(ll, k)) < 0.01


def test_penalty_ordering():
    ll, k, n = -100.0, 3, 100
    assert aic(ll, k) < hqic(ll, k, n) < bic(ll, k, n)


def test_more_parameters_penalized():
    ll, n = -100.0, 100
    assert aic(ll, 6) > aic(ll, 3)
    assert bic(ll, 6, n) > bic(ll, 3, n)


def test_gaussian_log_likelihood_matches_direct():
    n, rss = 50, 12.5
    sig2 = rss / n
    direct = -0.5 * n * (math.log(2 * math.pi) + math.log(sig2) + 1)
    assert abs(gaussian_log_likelihood(rss, n) - direct) < 1e-12


def test_better_fit_lower_aic():
    # Smaller RSS -> higher log-likelihood -> lower AIC at the same k.
    n, k = 50, 3
    good = aic(gaussian_log_likelihood(5.0, n), k)
    bad = aic(gaussian_log_likelihood(20.0, n), k)
    assert good < bad


def test_validation():
    with pytest.raises(ValueError):
        gaussian_log_likelihood(0.0, 10)
    with pytest.raises(ValueError):
        aicc(-100.0, 5, 5)          # n <= k + 1
    with pytest.raises(ValueError):
        hqic(-100.0, 3, 2)          # n < 3
    with pytest.raises(ValueError):
        aic(-100.0, -1)
