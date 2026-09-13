"""Two-dimensional numerical integration."""

import math

import pytest

from quantforge import integrate2d_gauss, integrate2d_simpson


def test_constant():
    assert abs(integrate2d_gauss(lambda x, y: 1.0, 0, 2, 0, 3) - 6.0) < 1e-9


def test_separable_polynomial():
    assert abs(integrate2d_gauss(lambda x, y: x * y, 0, 1, 0, 1) - 0.25) < 1e-12


def test_exact_for_polynomial_degree():
    # x^3 y^2 over [0,1]^2 = 1/4 * 1/3 = 1/12, exact at n=3 (degree 2n-1=5).
    assert abs(integrate2d_gauss(lambda x, y: x ** 3 * y ** 2, 0, 1, 0, 1, n=3)
               - 1.0 / 12.0) < 1e-12


def test_gaussian_via_simpson():
    r = integrate2d_simpson(lambda x, y: math.exp(-(x * x + y * y)),
                            -6, 6, -6, 6, nx=200, ny=200)
    assert abs(r - math.pi) < 1e-3


def test_sin_cos():
    r = integrate2d_simpson(lambda x, y: math.sin(x) * math.cos(y),
                            0, math.pi, 0, math.pi / 2, nx=100, ny=100)
    assert abs(r - 2.0) < 1e-4


def test_gauss_matches_simpson_on_polynomial():
    f = lambda x, y: x * x + y * y + x * y
    g = integrate2d_gauss(f, 0, 1, 0, 2, n=4)
    s = integrate2d_simpson(f, 0, 1, 0, 2, nx=100, ny=100)
    assert abs(g - s) < 1e-9


def test_plane():
    assert abs(integrate2d_gauss(lambda x, y: x + y, 0, 1, 0, 1, n=2) - 1.0) < 1e-12


def test_validation():
    with pytest.raises(ValueError):
        integrate2d_gauss(lambda x, y: 1.0, 0, 1, 0, 1, n=7)
