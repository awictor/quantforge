"""Tests for B-spline basis and curves, cross-checked against known properties and Bezier."""

import random

import pytest

from quantforge.bspline import bspline_basis, bspline_point, bspline_curve, open_uniform_knots
from quantforge.bezier import bezier_point


def _close(a, b, tol=1e-9):
    if isinstance(a, (int, float)):
        return abs(a - b) <= tol * max(1, abs(a), abs(b))
    return all(abs(x - y) <= tol * max(1, abs(x), abs(y)) for x, y in zip(a, b))


def test_fuzz_partition_of_unity_and_endpoints():
    rng = random.Random(371)
    for _ in range(3000):
        p = rng.randint(1, 4)
        n = rng.randint(p + 1, p + 6)
        knots = open_uniform_knots(n, p)
        for _ in range(4):
            t = rng.uniform(0, 0.999)
            s = sum(bspline_basis(i, p, knots, t) for i in range(n))
            assert _close(s, 1.0)
        ctrl = [rng.uniform(-10, 10) for _ in range(n)]
        assert _close(bspline_point(ctrl, p, 0.0), ctrl[0])
        assert _close(bspline_point(ctrl, p, 1.0), ctrl[-1])


def test_degree1_is_piecewise_linear_through_points():
    rng = random.Random(372)
    for _ in range(1000):
        n = rng.randint(2, 6)
        ctrl = [rng.uniform(-5, 5) for _ in range(n)]
        knots = open_uniform_knots(n, 1)
        for i in range(n):
            t = knots[i + 1]
            assert _close(bspline_point(ctrl, 1, t, knots), ctrl[i])


def test_bezier_equivalence_single_span():
    rng = random.Random(373)
    for _ in range(1000):
        n = rng.randint(2, 5)
        p = n - 1
        ctrl = [rng.uniform(-5, 5) for _ in range(n)]
        knots = open_uniform_knots(n, p)  # no interior knots
        for _ in range(3):
            t = rng.random()
            assert _close(bspline_point(ctrl, p, t, knots), bezier_point(ctrl, t), 1e-7)


def test_2d_endpoints_and_curve():
    c = [(0, 0), (1, 2), (3, 3), (4, 0)]
    assert _close(bspline_point(c, 3, 0.0), (0, 0))
    assert _close(bspline_point(c, 3, 1.0), (4, 0))
    curve = bspline_curve(c, 3, 10)
    assert _close(curve[0], (0, 0))
    assert _close(curve[-1], (4, 0))
    assert len(curve) == 10


def test_knot_vector_length_and_clamping():
    knots = open_uniform_knots(6, 3)
    assert len(knots) == 6 + 3 + 1
    assert knots[:4] == [0.0, 0.0, 0.0, 0.0]
    assert knots[-4:] == [1.0, 1.0, 1.0, 1.0]


def test_basis_nonnegative():
    knots = open_uniform_knots(5, 2)
    for i in range(5):
        for t in (0.0, 0.25, 0.5, 0.75, 0.999):
            assert bspline_basis(i, 2, knots, t) >= 0.0


def test_local_support():
    # moving one control point changes the curve only near that point's span
    n, p = 7, 3
    ctrl_a = [0.0] * n
    ctrl_b = list(ctrl_a)
    ctrl_b[n - 1] = 10.0  # bump the last point (support near t=1)
    # near t=0 the curve is unaffected by the last control point's basis
    assert _close(bspline_point(ctrl_a, p, 0.02), bspline_point(ctrl_b, p, 0.02))


def test_too_few_points_raises():
    with pytest.raises(ValueError):
        bspline_point([1, 2], 3, 0.5)


def test_knots_too_few_raises():
    with pytest.raises(ValueError):
        open_uniform_knots(2, 3)


def test_too_few_samples_raises():
    with pytest.raises(ValueError):
        bspline_curve([1, 2, 3], 1, 1)


def test_custom_knot_vector():
    ctrl = [0.0, 1.0, 2.0, 3.0]
    knots = open_uniform_knots(4, 2)
    v = bspline_point(ctrl, 2, 0.5, knots)
    assert isinstance(v, float)
