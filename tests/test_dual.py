"""Forward-mode autodiff with dual numbers."""

import cmath
import math

import pytest

from quantforge import Dual, dual_derivative
from quantforge.dual import exp, log, sqrt, sin, cos, tan, tanh
from quantforge.complex_step import complex_step_derivative


def test_elementary_derivatives_exact():
    cases = [
        (lambda x: x * x * x, lambda x: 3 * x ** 2, 2.0),
        (lambda x: exp(x), math.exp, 1.3),
        (lambda x: log(x), lambda x: 1 / x, 3.0),
        (lambda x: sin(x), math.cos, 0.7),
        (lambda x: sqrt(x), lambda x: 0.5 / math.sqrt(x), 4.0),
        (lambda x: x / (x + 1), lambda x: 1 / (x + 1) ** 2, 2.0),
        (lambda x: tanh(x), lambda x: 1 - math.tanh(x) ** 2, 0.5),
    ]
    for f, df, x in cases:
        assert abs(dual_derivative(f, x) - df(x)) < 1e-12


def test_chain_and_product_rule():
    f = lambda x: exp(x * x) * sin(x)
    df = lambda x: math.exp(x * x) * (2 * x * math.sin(x) + math.cos(x))
    assert abs(dual_derivative(f, 0.8) - df(0.8)) < 1e-12


def test_agrees_with_complex_step():
    d_dual = dual_derivative(lambda x: exp(x) + x * x * x, 1.1)
    d_cs = complex_step_derivative(lambda z: cmath.exp(z) + z ** 3, 1.1)
    assert abs(d_dual - d_cs) < 1e-12


def test_power_with_dual_exponent():
    f = lambda x: x ** x
    df = lambda x: x ** x * (math.log(x) + 1)
    assert abs(dual_derivative(f, 2.0) - df(2.0)) < 1e-10


def test_dual_arithmetic_value_and_deriv():
    x = Dual(3.0, 1.0)
    r = x * x + 2 * x + 1          # (x+1)^2, deriv 2(x+1) = 8 at x=3
    assert r.value == 16.0
    assert r.deriv == 8.0
