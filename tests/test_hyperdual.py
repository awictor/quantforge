"""Hyperdual numbers: exact first and second derivatives."""

import math

import pytest

from quantforge import HyperDual, hyperdual_derivatives, second_derivative
from quantforge.hyperdual import exp, log, sqrt, sin, cos
from quantforge.richardson_derivative import ridders_second_derivative


def test_first_and_second_derivatives_exact():
    cases = [
        (lambda x: x * x * x, lambda x: 3 * x ** 2, lambda x: 6 * x, 2.0),
        (lambda x: exp(x), math.exp, math.exp, 1.3),
        (lambda x: log(x), lambda x: 1 / x, lambda x: -1 / x ** 2, 3.0),
        (lambda x: sin(x), math.cos, lambda x: -math.sin(x), 0.7),
        (lambda x: sqrt(x), lambda x: 0.5 / math.sqrt(x), lambda x: -0.25 / x ** 1.5, 4.0),
        (lambda x: 1 / x, lambda x: -1 / x ** 2, lambda x: 2 / x ** 3, 2.5),
    ]
    for f, d1, d2, x in cases:
        v, fp, fpp = hyperdual_derivatives(f, x)
        assert abs(fp - d1(x)) < 1e-10
        assert abs(fpp - d2(x)) < 1e-10


def test_composite_second_derivative():
    f = lambda x: exp(x) * sin(x)
    _, _, fpp = hyperdual_derivatives(f, 0.8)
    assert abs(fpp - 2 * math.exp(0.8) * math.cos(0.8)) < 1e-12


def test_agrees_with_ridders_second():
    hd = second_derivative(lambda x: exp(x) + x ** 4, 1.1)
    rid, _ = ridders_second_derivative(lambda x: math.exp(x) + x ** 4, 1.1)
    assert abs(hd - rid) < 1e-6


def test_quotient_second_derivative():
    _, _, fpp = hyperdual_derivatives(lambda x: x / (x + 1), 2.0)
    assert abs(fpp - (-2 / 3 ** 3)) < 1e-12


def test_constant_second_derivative_is_two_for_square():
    assert abs(second_derivative(lambda x: x * x, 5.0) - 2.0) < 1e-12


def test_hyperdual_arithmetic():
    x = HyperDual(3.0, 1.0, 1.0, 0.0)
    r = x * x                         # value 9, f'=2x=6, f''=2
    assert r.f0 == 9.0
    assert r.f1 == 6.0
    assert r.f12 == 2.0
