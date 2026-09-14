"""Tests for Bezier curves / Bernstein polynomials, cross-checked against the Bernstein sum."""

import random

import pytest

from quantforge.bezier import (
    bernstein,
    bezier_point,
    bezier_curve,
    bezier_derivative_control,
    bezier_tangent,
    bezier_subdivide,
)


def _bernstein_eval(ctrl, t):
    n = len(ctrl) - 1
    if isinstance(ctrl[0], (int, float)):
        return sum(bernstein(n, i, t) * ctrl[i] for i in range(n + 1))
    dim = len(ctrl[0])
    return tuple(sum(bernstein(n, i, t) * ctrl[i][k] for i in range(n + 1)) for k in range(dim))


def _close(a, b, tol=1e-9):
    if isinstance(a, (int, float)):
        return abs(a - b) <= tol * max(1, abs(a), abs(b))
    return all(abs(x - y) <= tol * max(1, abs(x), abs(y)) for x, y in zip(a, b))


def test_fuzz_vs_bernstein_sum():
    rng = random.Random(351)
    for _ in range(4000):
        n = rng.randint(1, 6)
        dim = rng.choice([0, 2, 3])
        if dim == 0:
            ctrl = [rng.uniform(-10, 10) for _ in range(n + 1)]
        else:
            ctrl = [tuple(rng.uniform(-10, 10) for _ in range(dim)) for _ in range(n + 1)]
        for _ in range(4):
            t = rng.random()
            assert _close(bezier_point(ctrl, t), _bernstein_eval(ctrl, t))


def test_endpoints():
    ctrl = [(0, 0), (1, 5), (3, 2), (4, 4)]
    assert _close(bezier_point(ctrl, 0.0), ctrl[0])
    assert _close(bezier_point(ctrl, 1.0), ctrl[-1])


def test_derivative_vs_finite_difference():
    rng = random.Random(352)
    ctrl = [rng.uniform(-5, 5) for _ in range(5)]
    t = 0.4
    h = 1e-6
    fd = (bezier_point(ctrl, t + h) - bezier_point(ctrl, t - h)) / (2 * h)
    assert _close(bezier_tangent(ctrl, t), fd, 1e-4)


def test_subdivision_continuity_and_reparam():
    ctrl = [(0, 0), (1, 2), (3, 3), (5, 0)]
    ts = 0.37
    left, right = bezier_subdivide(ctrl, ts)
    assert _close(left[-1], bezier_point(ctrl, ts))
    assert _close(right[0], bezier_point(ctrl, ts))
    assert _close(left[0], ctrl[0])
    assert _close(right[-1], ctrl[-1])
    # a point on the left sub-curve at s equals the original at s*ts
    assert _close(bezier_point(left, 0.5), bezier_point(ctrl, 0.5 * ts), 1e-8)


def test_quadratic_explicit():
    assert _close(bezier_point([0, 1, 0], 0.5), 0.5)


def test_linear_is_lerp():
    assert _close(bezier_point([2, 8], 0.25), 3.5)


def test_2d_point():
    assert _close(bezier_point([(0, 0), (1, 2), (2, 0)], 0.5), (1.0, 1.0))


def test_bernstein_partition_of_unity():
    rng = random.Random(353)
    for n in range(6):
        for _ in range(5):
            t = rng.random()
            assert abs(sum(bernstein(n, i, t) for i in range(n + 1)) - 1.0) < 1e-12


def test_bernstein_out_of_range_is_zero():
    assert bernstein(3, -1, 0.5) == 0.0
    assert bernstein(3, 4, 0.5) == 0.0


def test_derivative_control_degree():
    ctrl = [0, 1, 4, 9]  # degree 3 -> derivative degree 2 (3 control points)
    d = bezier_derivative_control(ctrl)
    assert len(d) == 3
    assert d == [3, 9, 15]  # n*(P_{i+1}-P_i) = 3*(1,3,5)


def test_single_point_curve():
    assert bezier_point([7], 0.5) == 7
    assert bezier_tangent([7], 0.5) == 0.0
    assert bezier_derivative_control([7]) == []


def test_bezier_curve_samples():
    pts = bezier_curve([0, 10], 5)
    assert pts == [0.0, 2.5, 5.0, 7.5, 10.0]


def test_empty_control_raises():
    with pytest.raises(ValueError):
        bezier_point([], 0.5)


def test_too_few_samples_raises():
    with pytest.raises(ValueError):
        bezier_curve([1, 2], 1)
