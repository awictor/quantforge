"""Underwater curve, drawdown duration, and rolling Sharpe (perfmetrics)."""

import math

import pytest

from quantforge import (
    drawdown_curve, longest_drawdown_duration, rolling_sharpe,
    max_drawdown, sharpe_ratio,
)


R = [0.01, -0.005, 0.02, 0.008, -0.012, 0.015, 0.003, -0.002]


def test_curve_length_and_non_negative():
    curve = drawdown_curve(R)
    assert len(curve) == len(R)
    assert all(x >= 0.0 for x in curve)


def test_curve_max_equals_max_drawdown():
    assert max(drawdown_curve(R)) == pytest.approx(max_drawdown(R), abs=1e-12)


def test_curve_known_values():
    # equity 1.1, 0.55, 0.66 -> drawdowns 0, 0.5, 0.4.
    curve = drawdown_curve([0.1, -0.5, 0.2])
    assert curve == pytest.approx([0.0, 0.5, 0.4], abs=1e-12)


def test_longest_drawdown_duration():
    assert longest_drawdown_duration([0.1, -0.5, 0.2]) == 2  # never recovers
    assert longest_drawdown_duration([-0.1, 0.2]) == 1        # recovers period 2
    assert longest_drawdown_duration([0.01, 0.02]) == 0       # never underwater


def test_rolling_sharpe_length_and_first_value():
    rs = rolling_sharpe(R, 4)
    assert len(rs) == len(R) - 4 + 1
    assert rs[0] == pytest.approx(sharpe_ratio(R[:4]), abs=1e-12)


def test_rolling_sharpe_flat_window_is_nan():
    rs = rolling_sharpe([0.01, 0.01, 0.01, 0.02], 3)
    assert math.isnan(rs[0])   # first window is constant


def test_validation():
    with pytest.raises(ValueError):
        rolling_sharpe(R, 1)
    with pytest.raises(ValueError):
        rolling_sharpe(R, 100)
