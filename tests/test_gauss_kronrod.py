"""Adaptive Gauss-Kronrod quadrature."""

import math

import pytest

from quantforge import gauss_kronrod
from quantforge.gauss_kronrod import _gk15


def test_sine_and_exp():
    assert abs(gauss_kronrod(math.sin, 0.0, math.pi) - 2.0) < 1e-11
    assert abs(gauss_kronrod(math.exp, 0.0, 1.0) - (math.e - 1.0)) < 1e-11


def test_arctan_derivative():
    assert abs(gauss_kronrod(lambda x: 1 / (1 + x * x), 0.0, 1.0) - math.pi / 4) < 1e-11


def test_gaussian_matches_erf():
    val = gauss_kronrod(lambda x: math.exp(-x * x), 0.0, 2.0)
    assert abs(val - math.sqrt(math.pi) / 2 * math.erf(2)) < 1e-11


def test_sharp_peak_adaptive():
    peak = lambda x: 1.0 / ((x - 0.3) ** 2 + 0.001)
    exact = (math.atan((1 - 0.3) / 0.001 ** 0.5)
             - math.atan((0 - 0.3) / 0.001 ** 0.5)) / 0.001 ** 0.5
    assert abs(gauss_kronrod(peak, 0.0, 1.0, tol=1e-9) - exact) < 1e-4


def test_error_estimate_bounds_true_error():
    r, e = _gk15(math.exp, 0.0, 1.0)
    assert e >= abs(r - (math.e - 1.0)) - 1e-16


def test_reversed_limits_flip_sign():
    assert abs(gauss_kronrod(math.exp, 0.0, 1.0) + gauss_kronrod(math.exp, 1.0, 0.0)) < 1e-11


def test_polynomial_exact():
    assert abs(gauss_kronrod(lambda x: x ** 7, 0.0, 1.0) - 0.125) < 1e-12
