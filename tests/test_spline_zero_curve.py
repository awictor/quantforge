"""Cubic-spline-interpolated zero-rate discount curve."""

import math

import pytest

from quantforge import SplineZeroCurve


TIMES = [0.5, 1, 2, 5, 10]
ZEROS = [0.02, 0.025, 0.03, 0.035, 0.04]


def test_reprices_pillars():
    c = SplineZeroCurve(TIMES, ZEROS)
    assert all(c.zero_rate(TIMES[i]) == pytest.approx(ZEROS[i], abs=1e-9)
               for i in range(5))


def test_df_from_zero_rate():
    c = SplineZeroCurve(TIMES, ZEROS)
    assert c.df(2) == pytest.approx(math.exp(-0.03 * 2))
    assert c.df(0) == 1.0


def test_df_monotone_decreasing_upward_curve():
    c = SplineZeroCurve(TIMES, ZEROS)
    Ts = [0.5, 1, 2, 3, 5, 7, 10]
    dfs = [c.df(T) for T in Ts]
    assert all(dfs[i] > dfs[i + 1] for i in range(len(dfs) - 1))


def test_smooth_interpolation_between_pillars():
    c = SplineZeroCurve(TIMES, ZEROS)
    z = c.zero_rate(1.5)
    assert 0.025 < z < 0.03


def test_flat_extrapolation():
    c = SplineZeroCurve(TIMES, ZEROS)
    assert c.zero_rate(0.1) == ZEROS[0]
    assert c.zero_rate(20) == ZEROS[-1]


def test_forward_rate_positive():
    c = SplineZeroCurve(TIMES, ZEROS)
    assert c.forward_rate(1, 2) > 0


def test_validation():
    with pytest.raises(ValueError):
        SplineZeroCurve([1], [0.02])
    with pytest.raises(ValueError):
        SplineZeroCurve(TIMES, ZEROS).forward_rate(2, 1)
