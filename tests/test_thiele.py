"""Thiele's continued-fraction rational interpolation."""

import math

import pytest

from quantforge import thiele_coefficients, thiele_eval, thiele_interpolate


def test_passes_through_nodes():
    xs = [1, 2, 3, 4, 5]
    ys = [1 / (1 + x) for x in xs]
    for i in range(len(xs)):
        assert abs(thiele_interpolate(xs, ys, xs[i]) - ys[i]) < 1e-9


def test_recovers_rational_function():
    f = lambda x: (2 * x + 1) / (x * x + 1)
    xs = [0, 1, 2, 3, 4, -1]
    ys = [f(x) for x in xs]
    for xt in (0.5, 1.5, 2.7, -0.5, 5.0):
        assert abs(thiele_interpolate(xs, ys, xt) - f(xt)) < 1e-9


def test_tangent_with_poles():
    xs = [0.1 * i for i in range(1, 10)]
    ys = [math.tan(x) for x in xs]
    assert abs(thiele_interpolate(xs, ys, 0.55) - math.tan(0.55)) < 1e-6


def test_coefficients_then_eval():
    xs = [0, 1, 2, 3]
    ys = [1.0, 0.5, 1.0 / 3, 0.25]           # 1/(1+x)
    coef = thiele_coefficients(xs, ys)
    assert abs(thiele_eval(xs, coef, 1.5) - 1 / 2.5) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        thiele_coefficients([1], [1])
    with pytest.raises(ValueError):
        thiele_coefficients([1, 2], [3, 3])      # equal y at distinct x
