"""Gamma, chi-square, Poisson, F and binomial distributions."""

import math

import pytest

from quantforge import (gamma_cdf, gamma_pdf, gamma_ppf, chi2_cdf, chi2_sf,
                        chi2_ppf, poisson_pmf, poisson_cdf, f_cdf, f_ppf,
                        binomial_cdf, binomial_pmf)


def test_gamma_shape_one_is_exponential():
    for x in (0.5, 2.0, 5.0):
        assert abs(gamma_cdf(x, 1.0, 1.0) - (1.0 - math.exp(-x))) < 1e-12


def test_gamma_ppf_roundtrip():
    for p in (0.1, 0.5, 0.9):
        x = gamma_ppf(p, 2.5, 3.0)
        assert abs(gamma_cdf(x, 2.5, 3.0) - p) < 1e-9


def test_gamma_pdf_integrates_to_cdf():
    # Trapezoidal integral of the pdf matches the cdf.
    shape, scale, b = 3.0, 1.5, 8.0
    n = 20000
    h = b / n
    acc = 0.5 * (gamma_pdf(0.0, shape, scale) + gamma_pdf(b, shape, scale))
    for i in range(1, n):
        acc += gamma_pdf(i * h, shape, scale)
    assert abs(acc * h - gamma_cdf(b, shape, scale)) < 1e-4


def test_chi2_textbook_quantiles():
    assert abs(chi2_ppf(0.95, 1) - 3.841459) < 1e-4
    assert abs(chi2_ppf(0.95, 10) - 18.307038) < 1e-4
    assert abs(chi2_ppf(0.99, 5) - 15.086272) < 1e-4


def test_chi2_cdf_plus_sf_is_one():
    for x, df in [(5.0, 3), (1.0, 1), (20.0, 10)]:
        assert abs(chi2_cdf(x, df) + chi2_sf(x, df) - 1.0) < 1e-12


def test_poisson_cdf_matches_direct_sum():
    lam = 4.3
    for k in (0, 3, 7, 12):
        direct = sum(poisson_pmf(i, lam) for i in range(k + 1))
        assert abs(poisson_cdf(k, lam) - direct) < 1e-10


def test_poisson_pmf_sums_to_one():
    assert abs(sum(poisson_pmf(i, 4.3) for i in range(60)) - 1.0) < 1e-10


def test_poisson_zero_lambda():
    assert poisson_pmf(0, 0.0) == 1.0
    assert poisson_pmf(3, 0.0) == 0.0
    assert poisson_cdf(0, 0.0) == 1.0


def test_f_textbook_quantiles():
    assert abs(f_ppf(0.95, 1, 10) - 4.9646) < 1e-3
    assert abs(f_ppf(0.95, 5, 20) - 2.7109) < 1e-3


def test_f_ppf_roundtrip():
    for p in (0.25, 0.75, 0.9):
        x = f_ppf(p, 5, 10)
        assert abs(f_cdf(x, 5, 10) - p) < 1e-8


def test_binomial_cdf_matches_direct_sum():
    n, p = 20, 0.3
    for k in (3, 6, 10, 15):
        direct = sum(binomial_pmf(i, n, p) for i in range(k + 1))
        assert abs(binomial_cdf(k, n, p) - direct) < 1e-10


def test_binomial_pmf_sums_to_one():
    assert abs(sum(binomial_pmf(i, 20, 0.3) for i in range(21)) - 1.0) < 1e-10


def test_binomial_edge_probabilities():
    assert binomial_pmf(0, 10, 0.0) == 1.0
    assert binomial_pmf(10, 10, 1.0) == 1.0
    assert binomial_cdf(10, 10, 0.3) == 1.0


def test_validation():
    with pytest.raises(ValueError):
        gamma_cdf(1.0, 0.0)
    with pytest.raises(ValueError):
        gamma_ppf(1.5, 2.0)
    with pytest.raises(ValueError):
        chi2_cdf(1.0, 0.0)
    with pytest.raises(ValueError):
        poisson_pmf(1, -1.0)
    with pytest.raises(ValueError):
        f_cdf(1.0, 0.0, 5.0)
    with pytest.raises(ValueError):
        binomial_cdf(3, 10, 1.5)
