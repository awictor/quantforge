"""Chebyshev polynomial approximation."""

import math

import pytest

from quantforge import chebyshev_fit, chebyshev_eval, chebyshev_derivative


def test_exact_for_polynomial():
    poly = lambda x: 2 - 3 * x + x ** 2 - 0.5 * x ** 3
    c = chebyshev_fit(poly, -2, 3, 3)
    for x in (-1.5, 0.0, 1.0, 2.5):
        assert abs(chebyshev_eval(c, -2, 3, x) - poly(x)) < 1e-10


def test_spectral_accuracy_exp():
    c = chebyshev_fit(math.exp, 0, 2, 15)
    for x in (0.1, 0.5, 1.0, 1.9):
        assert abs(chebyshev_eval(c, 0, 2, x) - math.exp(x)) < 1e-12


def test_spectral_accuracy_sin():
    c = chebyshev_fit(math.sin, 0, math.pi, 20)
    for x in (0.3, 1.0, 2.0, 3.0):
        assert abs(chebyshev_eval(c, 0, math.pi, x) - math.sin(x)) < 1e-12


def test_derivative_of_exp():
    c = chebyshev_fit(math.exp, 0, 2, 15)
    dc = chebyshev_derivative(c, 0, 2)
    for x in (0.1, 0.5, 1.0, 1.9):
        assert abs(chebyshev_eval(dc, 0, 2, x) - math.exp(x)) < 1e-10


def test_derivative_of_sin_is_cos():
    c = chebyshev_fit(math.sin, 0, math.pi, 20)
    dc = chebyshev_derivative(c, 0, math.pi)
    for x in (0.3, 1.0, 2.0, 3.0):
        assert abs(chebyshev_eval(dc, 0, math.pi, x) - math.cos(x)) < 1e-10


def test_coefficients_decay_for_smooth():
    c = chebyshev_fit(math.exp, 0, 2, 20)
    # Later coefficients are much smaller than the leading ones.
    assert abs(c[-1]) < 1e-12
    assert abs(c[0]) > abs(c[10])


def test_validation():
    with pytest.raises(ValueError):
        chebyshev_fit(math.exp, 0, 1, 0)
    with pytest.raises(ValueError):
        chebyshev_eval([], 0, 1, 0.5)
