"""General Richardson extrapolation."""

import math

import pytest

from quantforge import richardson_extrapolate, richardson_table


def test_forward_difference_first_order():
    f = lambda h: (math.exp(h) - 1) / h
    est = [f(0.1 / 2 ** i) for i in range(6)]
    assert abs(richardson_extrapolate(est, p=1) - 1.0) < 1e-9


def test_central_difference_second_order():
    cd = lambda h: (math.sin(1 + h) - math.sin(1 - h)) / (2 * h)
    est = [cd(0.2 / 2 ** i) for i in range(6)]
    assert abs(richardson_extrapolate(est, p=2) - math.cos(1)) < 1e-10


def test_romberg_from_trapezoid():
    def trap(n):
        a, b = 0.0, math.pi
        h = (b - a) / n
        return h * (0.5 * math.sin(a) + 0.5 * math.sin(b)
                    + sum(math.sin(a + i * h) for i in range(1, n)))
    est = [trap(2 ** i) for i in range(1, 7)]
    assert abs(richardson_extrapolate(est, p=2) - 2.0) < 1e-9


def test_exact_on_error_model():
    est = [5 + 2 * (0.1 / 2 ** i) for i in range(4)]
    assert abs(richardson_extrapolate(est, p=1) - 5.0) < 1e-12


def test_table_diagonal_converges():
    f = lambda h: (math.exp(h) - 1) / h
    est = [f(0.4 / 2 ** i) for i in range(5)]
    T = richardson_table(est, p=1)
    assert abs(T[-1][-1] - 1.0) < 1e-6


def test_single_estimate():
    assert richardson_extrapolate([3.14]) == 3.14


def test_validation():
    with pytest.raises(ValueError):
        richardson_extrapolate([])
    with pytest.raises(ValueError):
        richardson_table([])
