"""Tests for the Adams-Bashforth-Moulton ODE integrator, vs closed-form solutions."""

import math

import pytest

from quantforge.multistep_ode import adams_bashforth_moulton


def _close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def test_exponential_growth():
    ts, ys = adams_bashforth_moulton(lambda t, y: y, 1.0, 0, 1, 200)
    assert _close(ys[-1], math.e, 1e-6)


def test_exponential_decay():
    ts, ys = adams_bashforth_moulton(lambda t, y: -2 * y, 3.0, 0, 2, 400)
    assert _close(ys[-1], 3 * math.exp(-4), 1e-6)


def test_linear_rhs():
    # y' = t, y(0)=0 -> t^2/2; at t=3 -> 4.5
    ts, ys = adams_bashforth_moulton(lambda t, y: t, 0.0, 0, 3, 300)
    assert _close(ys[-1], 4.5, 1e-8)


def test_harmonic_system():
    # y0'=y1, y1'=-y0 from (0,1) -> (sin t, cos t)
    def harm(t, y):
        return [y[1], -y[0]]
    ts, ys = adams_bashforth_moulton(harm, [0.0, 1.0], 0, math.pi, 500)
    assert _close(ys[-1][0], math.sin(math.pi), 1e-5)
    assert _close(ys[-1][1], math.cos(math.pi), 1e-5)


def test_fourth_order_convergence():
    def err(n):
        _, ys = adams_bashforth_moulton(lambda t, y: y, 1.0, 0, 1, n)
        return abs(ys[-1] - math.e)
    e1 = err(50)
    e2 = err(100)
    assert e1 / e2 > 8  # ~16 for a 4th-order method


def test_time_points():
    ts, ys = adams_bashforth_moulton(lambda t, y: y, 1.0, 0.0, 2.0, 4)
    assert len(ts) == 5
    assert ts[0] == 0.0
    assert _close(ts[-1], 2.0)
    assert len(ys) == 5


def test_short_integration_uses_rk4():
    # fewer than 4 steps: bootstrap only, still accurate
    ts, ys = adams_bashforth_moulton(lambda t, y: y, 1.0, 0, 1, 2)
    assert len(ys) == 3
    assert _close(ys[-1], math.e, 1e-3)


def test_vector_system_logistic_pair():
    # two decoupled decays
    def f(t, y):
        return [-y[0], -3 * y[1]]
    ts, ys = adams_bashforth_moulton(f, [2.0, 5.0], 0, 1, 300)
    assert _close(ys[-1][0], 2 * math.exp(-1), 1e-6)
    assert _close(ys[-1][1], 5 * math.exp(-3), 1e-6)


def test_constant_solution():
    # y' = 0 -> constant
    ts, ys = adams_bashforth_moulton(lambda t, y: 0.0, 7.0, 0, 5, 10)
    assert all(_close(y, 7.0) for y in ys)


def test_single_step():
    ts, ys = adams_bashforth_moulton(lambda t, y: y, 1.0, 0, 0.1, 1)
    assert len(ys) == 2
    assert _close(ys[-1], math.exp(0.1), 1e-4)


def test_n_steps_below_one_raises():
    with pytest.raises(ValueError):
        adams_bashforth_moulton(lambda t, y: y, 1.0, 0, 1, 0)
