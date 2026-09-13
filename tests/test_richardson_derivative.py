"""High-accuracy Ridders numerical derivatives."""

import math

import pytest

from quantforge import ridders_derivative, ridders_second_derivative, price, delta


def test_first_derivatives_near_machine_precision():
    cases = [
        (math.sin, math.cos, 0.7),
        (math.exp, math.exp, 1.3),
        (lambda x: x ** 4, lambda x: 4 * x ** 3, 2.0),
        (math.log, lambda x: 1 / x, 3.0),
        (lambda x: 1 / x, lambda x: -1 / x ** 2, 2.5),
    ]
    for f, df, x in cases:
        d, err = ridders_derivative(f, x)
        assert abs(d - df(x)) < 1e-9


def test_second_derivatives():
    cases = [
        (math.sin, lambda x: -math.sin(x), 0.7),
        (math.exp, math.exp, 1.3),
        (lambda x: x ** 4, lambda x: 12 * x ** 2, 2.0),
    ]
    for f, d2, x in cases:
        d, err = ridders_second_derivative(f, x)
        assert abs(d - d2(x)) < 1e-6


def test_error_estimate_is_honest():
    d, err = ridders_derivative(math.sin, 0.7)
    assert abs(d - math.cos(0.7)) <= 20 * err


def test_matches_black_scholes_delta():
    d, _ = ridders_derivative(lambda S: price(S, 100, 1.0, 0.05, 0.2), 100.0)
    assert abs(d - delta(100, 100, 1.0, 0.05, 0.2)) < 1e-6


def test_validation():
    with pytest.raises(ValueError):
        ridders_derivative(math.sin, 0.5, h=0)
    with pytest.raises(ValueError):
        ridders_second_derivative(math.sin, 0.5, h=-1)
