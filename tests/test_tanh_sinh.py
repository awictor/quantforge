"""Tanh-sinh (double-exponential) quadrature, including singular integrands."""

import math

import pytest

from quantforge import tanh_sinh


def test_inverse_sqrt_endpoint_singularity():
    # integral_0^1 x^{-1/2} dx = 2
    assert abs(tanh_sinh(lambda x: 1.0 / math.sqrt(x), 0.0, 1.0, levels=6) - 2.0) < 1e-7


def test_log_singularity():
    # integral_0^1 -ln(x) dx = 1
    assert abs(tanh_sinh(lambda x: -math.log(x), 0.0, 1.0, levels=6) - 1.0) < 1e-8


def test_log_over_sqrt():
    # integral_0^1 ln(x)/sqrt(x) dx = -4
    assert abs(tanh_sinh(lambda x: math.log(x) / math.sqrt(x), 0.0, 1.0, levels=7) + 4.0) < 1e-4


def test_semicircle_two_sided_singularity():
    # integral_-1^1 sqrt(1 - x^2) dx = pi/2
    assert abs(tanh_sinh(lambda x: math.sqrt(1 - x * x), -1.0, 1.0, levels=7)
               - math.pi / 2) < 1e-7


def test_smooth_polynomial_exact():
    assert abs(tanh_sinh(lambda x: x * x, 0.0, 1.0, levels=5) - 1.0 / 3.0) < 1e-9


def test_gaussian_matches_erf():
    val = tanh_sinh(lambda x: math.exp(-x * x), 0.0, 2.0, levels=6)
    assert abs(val - math.sqrt(math.pi) / 2 * math.erf(2)) < 1e-8


def test_reversed_limits_flip_sign():
    fwd = tanh_sinh(lambda x: math.exp(x), 0.0, 1.0, levels=5)
    rev = tanh_sinh(lambda x: math.exp(x), 1.0, 0.0, levels=5)
    assert abs(fwd + rev) < 1e-9


def test_more_levels_improve_accuracy():
    exact = 2.0
    coarse = abs(tanh_sinh(lambda x: 1.0 / math.sqrt(x), 0.0, 1.0, levels=3) - exact)
    fine = abs(tanh_sinh(lambda x: 1.0 / math.sqrt(x), 0.0, 1.0, levels=7) - exact)
    assert fine < coarse


def test_validation():
    with pytest.raises(ValueError):
        tanh_sinh(lambda x: x, 0.0, 1.0, levels=-1)
