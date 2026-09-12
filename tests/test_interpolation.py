"""Cubic interpolation: natural spline and monotone Hermite."""

import pytest

from quantforge import natural_cubic_spline, monotone_cubic


XS = [0, 1, 2, 3, 4]
YS = [0, 1, 4, 9, 16]


def test_spline_passes_through_knots():
    f = natural_cubic_spline(XS, YS)
    assert all(f(XS[i]) == pytest.approx(YS[i], abs=1e-9) for i in range(5))


def test_spline_exact_on_linear():
    f = natural_cubic_spline([0, 1, 2, 3], [0, 2, 4, 6])
    assert f(1.5) == pytest.approx(3.0)
    assert f(0.5) == pytest.approx(1.0)


def test_monotone_passes_through_knots():
    g = monotone_cubic(XS, YS)
    assert all(g(XS[i]) == pytest.approx(YS[i], abs=1e-9) for i in range(5))


def test_monotone_preserves_monotonicity():
    g = monotone_cubic(XS, YS)
    pts = [g(x / 10) for x in range(0, 41)]
    assert all(pts[i] <= pts[i + 1] + 1e-9 for i in range(len(pts) - 1))


def test_monotone_no_overshoot_where_spline_does():
    xm, ym = [0, 1, 2, 3], [0, 0, 0, 1]
    g = monotone_cubic(xm, ym)
    vals = [g(x / 10) for x in range(0, 31)]
    assert all(-1e-9 <= v <= 1 + 1e-9 for v in vals)
    # The natural spline overshoots on the same data.
    f = natural_cubic_spline(xm, ym)
    svals = [f(x / 10) for x in range(0, 31)]
    assert min(svals) < -1e-6 or max(svals) > 1 + 1e-6


def test_monotone_exact_on_linear():
    g = monotone_cubic([0, 1, 2, 3], [0, 2, 4, 6])
    assert g(1.5) == pytest.approx(3.0)


def test_validation():
    with pytest.raises(ValueError):
        natural_cubic_spline([0], [0])
    with pytest.raises(ValueError):
        monotone_cubic([0, 0, 1], [0, 1, 2])   # non-increasing xs
