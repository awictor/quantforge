"""Gauss-Hermite quadrature for Gaussian-weighted expectations."""

import math

import pytest

from quantforge import gauss_hermite_nodes_weights, gauss_hermite_expectation


def test_weights_sum_to_one():
    for n in (1, 2, 4, 8, 16, 24):
        _, w = gauss_hermite_nodes_weights(n)
        assert abs(sum(w) - 1.0) < 1e-12


def test_nodes_symmetric_and_sorted():
    x, _ = gauss_hermite_nodes_weights(8)
    assert all(x[i] < x[i + 1] for i in range(len(x) - 1))
    for i in range(len(x)):
        assert abs(x[i] + x[-1 - i]) < 1e-9


def test_standard_normal_moments():
    x, w = gauss_hermite_nodes_weights(8)
    m2 = sum(w[i] * x[i] ** 2 for i in range(8))
    m4 = sum(w[i] * x[i] ** 4 for i in range(8))
    m6 = sum(w[i] * x[i] ** 6 for i in range(8))
    assert abs(m2 - 1.0) < 1e-10
    assert abs(m4 - 3.0) < 1e-9
    assert abs(m6 - 15.0) < 1e-8


def test_odd_moments_vanish():
    x, w = gauss_hermite_nodes_weights(10)
    for p in (1, 3, 5):
        assert abs(sum(w[i] * x[i] ** p for i in range(10))) < 1e-9


def test_exact_for_polynomial_up_to_degree_2n_minus_1():
    # n=3 is exact through degree 5; E[Z^4]=3 (deg 4) must be exact.
    x, w = gauss_hermite_nodes_weights(3)
    assert abs(sum(w[i] * x[i] ** 4 for i in range(3)) - 3.0) < 1e-10


def test_mgf_matches_closed_form():
    for sig in (0.2, 0.5, 1.0):
        approx = gauss_hermite_expectation(lambda z: math.exp(sig * z), 0.0, 1.0, n=24)
        assert abs(approx - math.exp(sig * sig / 2)) < 1e-7


def test_expectation_shift_and_scale():
    # E[X^2] for X ~ N(mu, sigma^2) is mu^2 + sigma^2.
    mu, sigma = 1.5, 0.7
    val = gauss_hermite_expectation(lambda x: x * x, mu, sigma, n=8)
    assert abs(val - (mu * mu + sigma * sigma)) < 1e-9


def test_lognormal_mean():
    # E[exp(X)] for X ~ N(mu, sigma^2) is exp(mu + sigma^2/2).
    mu, sigma = 0.05, 0.3
    val = gauss_hermite_expectation(math.exp, mu, sigma, n=24)
    assert abs(val - math.exp(mu + sigma * sigma / 2)) < 1e-8


def test_zero_sigma_collapses_to_point():
    assert gauss_hermite_expectation(lambda x: x ** 3 + 2, mu=2.0, sigma=0.0) == 10.0


def test_validation():
    with pytest.raises(ValueError):
        gauss_hermite_nodes_weights(0)
    with pytest.raises(ValueError):
        gauss_hermite_expectation(lambda x: x, sigma=-1.0)
