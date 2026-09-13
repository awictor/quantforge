"""Clenshaw-Curtis quadrature."""

import math

import pytest

from quantforge import clenshaw_curtis


def test_sine_over_half_period():
    assert abs(clenshaw_curtis(math.sin, 0.0, math.pi) - 2.0) < 1e-12


def test_exponential():
    assert abs(clenshaw_curtis(math.exp, 0.0, 1.0) - (math.e - 1.0)) < 1e-12


def test_polynomial_exact_at_matching_order():
    # A degree-7 polynomial is integrated exactly by n = 8.
    assert abs(clenshaw_curtis(lambda x: x ** 7, 0.0, 1.0, n=8) - 0.125) < 1e-12


def test_arctan_derivative():
    assert abs(clenshaw_curtis(lambda x: 1 / (1 + x * x), 0.0, 1.0) - math.pi / 4) < 1e-12


def test_gaussian_matches_erf():
    val = clenshaw_curtis(lambda x: math.exp(-x * x), 0.0, 2.0)
    assert abs(val - math.sqrt(math.pi) / 2 * math.erf(2)) < 1e-12


def test_runge_function():
    # The Runge function is hard for equispaced rules; Chebyshev points handle it.
    val = clenshaw_curtis(lambda x: 1 / (1 + 25 * x * x), -1.0, 1.0, n=128)
    assert abs(val - 2.0 / 5.0 * math.atan(5)) < 1e-6


def test_reversed_limits_flip_sign():
    assert abs(clenshaw_curtis(math.exp, 0.0, 1.0) + clenshaw_curtis(math.exp, 1.0, 0.0)) < 1e-12


def test_agrees_with_romberg():
    from quantforge import romberg
    f = lambda x: math.cos(x) * math.exp(-0.1 * x)
    assert abs(clenshaw_curtis(f, 0.0, 3.0) - romberg(f, 0.0, 3.0)) < 1e-9


def test_odd_n_rounded_up():
    # n = 7 is rounded to 8; result matches an explicit even call.
    a = clenshaw_curtis(math.sin, 0.0, math.pi, n=7)
    b = clenshaw_curtis(math.sin, 0.0, math.pi, n=8)
    assert abs(a - b) < 1e-15


def test_validation():
    with pytest.raises(ValueError):
        clenshaw_curtis(math.sin, 0.0, 1.0, n=1)
