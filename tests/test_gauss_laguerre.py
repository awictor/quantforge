"""Gauss-Laguerre quadrature over the half-line."""

import math

import pytest

from quantforge import gauss_laguerre_nodes_weights, gauss_laguerre_integral


def test_weights_sum_to_one():
    for n in (1, 2, 4, 8, 16, 32):
        _, w = gauss_laguerre_nodes_weights(n)
        assert abs(sum(w) - 1.0) < 1e-12


def test_nodes_positive_and_sorted():
    x, _ = gauss_laguerre_nodes_weights(12)
    assert all(xi > 0.0 for xi in x)
    assert all(x[i] < x[i + 1] for i in range(len(x) - 1))


def test_moments_are_factorials():
    x, w = gauss_laguerre_nodes_weights(10)
    for m in range(6):
        approx = sum(w[i] * x[i] ** m for i in range(10))
        assert abs(approx - math.factorial(m)) < 5e-6


def test_exact_for_polynomial_up_to_degree_2n_minus_1():
    # n=3 exact through degree 5; the degree-4 moment 4! = 24 must be exact.
    x, w = gauss_laguerre_nodes_weights(3)
    assert abs(sum(w[i] * x[i] ** 4 for i in range(3)) - 24.0) < 1e-9


def test_exponential_decay_integral():
    # integral_0^inf e^{-2x} dx = 1/2
    val = gauss_laguerre_integral(lambda x: math.exp(-2 * x), n=32, rate=2.0)
    assert abs(val - 0.5) < 1e-9


def test_polynomial_times_exponential():
    # integral_0^inf x^2 e^{-3x} dx = 2 / 27
    val = gauss_laguerre_integral(lambda x: x * x * math.exp(-3 * x), n=40, rate=3.0)
    assert abs(val - 2.0 / 27.0) < 1e-8


def test_gaussian_half_line():
    # integral_0^inf e^{-x^2} dx = sqrt(pi) / 2
    val = gauss_laguerre_integral(lambda x: math.exp(-x * x), n=60, rate=1.0)
    assert abs(val - math.sqrt(math.pi) / 2) < 1e-5


def test_exponential_integral_form():
    # integral_0^inf e^{-x} / (1 + x) dx = 0.5963473623...
    val = gauss_laguerre_integral(lambda x: math.exp(-x) / (1 + x), n=48, rate=1.0)
    assert abs(val - 0.5963473623231941) < 1e-5


def test_validation():
    with pytest.raises(ValueError):
        gauss_laguerre_nodes_weights(0)
    with pytest.raises(ValueError):
        gauss_laguerre_integral(lambda x: math.exp(-x), rate=0.0)
