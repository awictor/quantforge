"""Romberg integration (Richardson extrapolation on the trapezoid rule)."""

import math

import pytest

from quantforge import romberg


def test_sine_over_half_period():
    assert abs(romberg(math.sin, 0.0, math.pi) - 2.0) < 1e-12


def test_exponential():
    assert abs(romberg(math.exp, 0.0, 1.0) - (math.e - 1.0)) < 1e-12


def test_polynomial_exact():
    assert abs(romberg(lambda x: x ** 7, 0.0, 1.0) - 0.125) < 1e-12


def test_arctan_derivative_gives_pi_over_4():
    assert abs(romberg(lambda x: 1.0 / (1.0 + x * x), 0.0, 1.0) - math.pi / 4) < 1e-12


def test_gaussian_matches_erf():
    val = romberg(lambda x: math.exp(-x * x), 0.0, 2.0)
    assert abs(val - math.sqrt(math.pi) / 2 * math.erf(2)) < 1e-12


def test_reversed_limits_flip_sign():
    assert abs(romberg(math.exp, 0.0, 1.0) + romberg(math.exp, 1.0, 0.0)) < 1e-12


def test_agrees_with_simpson_on_smooth():
    from quantforge import simpson
    f = lambda x: math.cos(x) * math.exp(-0.1 * x)
    assert abs(romberg(f, 0.0, 3.0) - simpson(f, 0.0, 3.0, 10000)) < 1e-8


def test_low_order_still_reasonable():
    # A single Richardson step already beats the raw trapezoid.
    approx = romberg(math.sin, 0.0, math.pi, max_order=2, tol=0.0)
    assert abs(approx - 2.0) < 1e-2


def test_validation():
    with pytest.raises(ValueError):
        romberg(math.sin, 0.0, 1.0, max_order=0)
