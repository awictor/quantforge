"""ODE solvers: RK4 and adaptive RK45."""

import math

import pytest

from quantforge import rk4, rk45


def test_rk4_exponential():
    ts, ys = rk4(lambda t, y: y, 0, 1.0, 1.0, n=100)
    assert abs(ys[-1][0] - math.e) < 1e-8


def test_rk45_exponential():
    ts, ys = rk45(lambda t, y: y, 0, 1.0, 1.0, tol=1e-10)
    assert abs(ys[-1][0] - math.e) < 1e-8


def test_decaying_exponential():
    ts, ys = rk4(lambda t, y: -y, 0, 1.0, 3.0, n=300)
    assert abs(ys[-1][0] - math.exp(-3)) < 1e-7


def test_harmonic_oscillator_system():
    ts, ys = rk45(lambda t, y: [y[1], -y[0]], 0, [1.0, 0.0], 2 * math.pi, tol=1e-10)
    assert abs(ys[-1][0] - 1.0) < 1e-5
    assert abs(ys[-1][1]) < 1e-5


def test_logistic_matches_closed_form():
    ts, ys = rk4(lambda t, y: y * (1 - y), 0, 0.1, 5.0, n=500)
    true = 1 / (1 + 9 * math.exp(-5))
    assert abs(ys[-1][0] - true) < 1e-5


def test_rk45_sin():
    ts, ys = rk45(lambda t, y: math.cos(t), 0, 0.0, 10.0, tol=1e-8)
    assert abs(ys[-1][0] - math.sin(10)) < 1e-6


def test_scalar_vector_consistency():
    _, ysc = rk4(lambda t, y: y, 0, 1.0, 1.0, n=50)
    _, yv = rk4(lambda t, y: [y[0]], 0, [1.0], 1.0, n=50)
    assert abs(ysc[-1][0] - yv[-1][0]) < 1e-15


def test_validation():
    with pytest.raises(ValueError):
        rk4(lambda t, y: y, 0, 1.0, 1.0, n=0)
