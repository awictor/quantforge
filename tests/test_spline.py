"""Tests for the natural cubic spline and smile interpolation."""

import math

import pytest

from quantforge import CubicSpline, SmileSpline


def test_interpolates_nodes_exactly():
    xs = [80, 90, 100, 110, 120]
    ys = [0.28, 0.24, 0.22, 0.225, 0.24]
    s = CubicSpline(xs, ys)
    for x, y in zip(xs, ys):
        assert s(x) == pytest.approx(y, abs=1e-12)


def test_reproduces_a_line():
    # A natural spline through collinear points is that line.
    s = CubicSpline([0, 1, 2, 3, 4], [1, 3, 5, 7, 9])
    for x in (0.5, 1.7, 2.3, 3.9):
        assert s(x) == pytest.approx(1 + 2 * x, abs=1e-9)


def test_reproduces_a_cubic_in_the_interior():
    # A smooth cubic is recovered well by the spline away from the ends, where
    # the natural (zero second-derivative) boundary intentionally deviates from
    # a genuine cubic's curvature.
    f = lambda x: x ** 3 - 2 * x + 1
    xs = [i * 0.5 for i in range(-6, 7)]  # -3 .. 3
    s = CubicSpline(xs, [f(x) for x in xs])
    for x in (-0.4, 0.9, 1.3):  # interior points
        assert s(x) == pytest.approx(f(x), abs=0.05)


def test_natural_boundary_second_derivative_zero():
    xs = [0, 1, 2, 3, 4]
    ys = [0, 1, 0, 1, 0]
    s = CubicSpline(xs, ys)
    assert s.m[0] == pytest.approx(0.0)
    assert s.m[-1] == pytest.approx(0.0)


def test_continuous_across_segments():
    # No jumps at the interior knots.
    xs = [0, 1, 2, 3]
    ys = [0, 2, -1, 3]
    s = CubicSpline(xs, ys)
    for k in (1, 2):
        left = s(k - 1e-7)
        right = s(k + 1e-7)
        assert left == pytest.approx(right, abs=1e-5)


def test_clamps_outside_range():
    s = CubicSpline([1, 2, 3], [10, 20, 30])
    assert s(0) == pytest.approx(10.0)   # below range -> first value
    assert s(5) == pytest.approx(30.0)   # above range -> last value


def test_smile_flat_extrapolation():
    strikes = [80, 90, 100, 110, 120]
    vols = [0.28, 0.24, 0.22, 0.225, 0.24]
    sm = SmileSpline(strikes, vols)
    assert sm.vol(80) == pytest.approx(0.28)
    assert sm.vol(50) == pytest.approx(0.28)   # flat below
    assert sm.vol(200) == pytest.approx(0.24)  # flat above
    # Interior value lies between neighboring quotes.
    v = sm.vol(95)
    assert 0.22 <= v <= 0.24


def test_requires_three_points():
    with pytest.raises(ValueError):
        CubicSpline([1, 2], [3, 4])


def test_rejects_non_increasing_x():
    with pytest.raises(ValueError):
        CubicSpline([1, 1, 2], [1, 2, 3])
