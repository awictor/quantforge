"""Numerical integration routines."""

import math

import pytest

from quantforge import trapezoid, simpson, gauss_legendre, adaptive_simpson


def test_trapezoid_polynomial():
    assert trapezoid(lambda x: x * x, 0, 1, 1000) == pytest.approx(1 / 3, abs=1e-5)


def test_simpson_exact_for_cubic():
    assert simpson(lambda x: x ** 3, 0, 2, 10) == pytest.approx(4.0, abs=1e-10)


def test_gauss_legendre_exact_to_degree_2n_minus_1():
    assert gauss_legendre(lambda x: x ** 5, 0, 1, 3) == pytest.approx(1 / 6, abs=1e-10)
    assert gauss_legendre(lambda x: x ** 9, 0, 1, 5) == pytest.approx(1 / 10, abs=1e-9)


def test_adaptive_simpson_sin():
    assert adaptive_simpson(math.sin, 0, math.pi) == pytest.approx(2.0, abs=1e-9)


def test_adaptive_simpson_peaked():
    assert adaptive_simpson(lambda x: 1 / (1 + x * x), -5, 5) == pytest.approx(
        2 * math.atan(5), abs=1e-8)


def test_adaptive_simpson_gaussian():
    exact = math.sqrt(math.pi) * math.erf(3)
    assert adaptive_simpson(lambda x: math.exp(-x * x), -3, 3) == pytest.approx(
        exact, abs=1e-8)


def test_validation():
    with pytest.raises(ValueError):
        gauss_legendre(math.sin, 0, 1, 7)
    with pytest.raises(ValueError):
        trapezoid(math.sin, 0, 1, 0)
