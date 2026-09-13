"""Complex-step differentiation."""

import cmath
import math

import pytest

from quantforge import complex_step_derivative, complex_step_gradient
from quantforge import ridders_derivative


def test_derivatives_machine_precision():
    cases = [
        (cmath.sin, math.cos, 0.7),
        (cmath.exp, math.exp, 1.3),
        (lambda z: z ** 4, lambda x: 4 * x ** 3, 2.0),
        (cmath.log, lambda x: 1 / x, 3.0),
        (lambda z: 1 / z, lambda x: -1 / x ** 2, 2.5),
        (lambda z: z ** 0.5, lambda x: 0.5 * x ** -0.5, 4.0),
    ]
    for f, df, x in cases:
        assert abs(complex_step_derivative(f, x) - df(x)) < 1e-12


def test_tiny_step_no_cancellation():
    # A finite difference at h=1e-100 would be pure round-off; complex-step is exact.
    d = complex_step_derivative(cmath.sin, 0.7, h=1e-100)
    assert abs(d - math.cos(0.7)) < 1e-12


def test_agrees_with_ridders():
    r, _ = ridders_derivative(math.sin, 0.7)
    c = complex_step_derivative(cmath.sin, 0.7)
    assert abs(r - c) < 1e-10


def test_gradient():
    g = lambda v: v[0] ** 2 + 3 * v[1] ** 2 + v[0] * v[1]
    grad = complex_step_gradient(g, [1.0, 2.0])
    assert abs(grad[0] - 4.0) < 1e-12       # 2x0 + x1
    assert abs(grad[1] - 13.0) < 1e-12      # 6x1 + x0


def test_validation():
    with pytest.raises(ValueError):
        complex_step_derivative(cmath.sin, 0.5, h=0)
