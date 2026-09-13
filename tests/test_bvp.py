"""Shooting method for two-point boundary-value problems."""

import math

import pytest

from quantforge import shooting_bvp


def test_sinh_problem():
    # y'' = y, y(0)=0, y(1)=1 -> y = sinh(t)/sinh(1), slope y'(0) = 1/sinh(1).
    r = shooting_bvp(lambda t, y, yp: y, 0, 1, 0.0, 1.0, -5, 5)
    assert abs(r["slope"] - 1 / math.sinh(1)) < 1e-6
    assert abs(r["ys"][-1][0] - 1.0) < 1e-6


def test_sine_problem():
    # y'' = -y, y(0)=0, y(pi/2)=1 -> y = sin(t), slope 1.
    r = shooting_bvp(lambda t, y, yp: -y, 0, math.pi / 2, 0.0, 1.0, -5, 5)
    assert abs(r["slope"] - 1.0) < 1e-6


def test_cubic_problem():
    # y'' = 6t, y(0)=0, y(1)=1 -> y = t^3, slope y'(0) = 0.
    r = shooting_bvp(lambda t, y, yp: 6 * t, 0, 1, 0.0, 1.0, -5, 5)
    assert abs(r["slope"]) < 1e-6
    # Every computed point lies on t^3.
    for t, y in zip(r["ts"], r["ys"]):
        assert abs(y[0] - t ** 3) < 1e-5


def test_terminal_boundary_hit():
    r = shooting_bvp(lambda t, y, yp: 6 * t, 0, 1, 0.0, 1.0, -5, 5)
    assert abs(r["ys"][-1][0] - 1.0) < 1e-7


def test_with_first_derivative_term():
    # y'' = -y', y(0)=0, y(1)=1.
    r = shooting_bvp(lambda t, y, yp: -yp, 0, 1, 0.0, 1.0, -5, 5)
    assert abs(r["ys"][-1][0] - 1.0) < 1e-6


def test_bad_bracket_raises():
    with pytest.raises(ValueError):
        # Both slopes overshoot the same way -> no sign change.
        shooting_bvp(lambda t, y, yp: y, 0, 1, 0.0, 1.0, 10, 20)
