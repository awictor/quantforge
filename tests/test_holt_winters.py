"""Exponential smoothing: Holt and Holt-Winters."""

import pytest

from quantforge import holt_linear, holt_winters_add


def test_holt_recovers_linear_trend():
    y = [float(i) for i in range(20)]      # slope 1
    level, trend, fc = holt_linear(y, 0.5, 0.3, 5)
    assert abs(trend - 1.0) < 0.05
    assert all(abs(fc[i] - (19 + i + 1)) < 0.5 for i in range(5))


def test_holt_forecast_is_straight_line():
    y = [float(i) for i in range(20)]
    _, _, fc = holt_linear(y, 0.5, 0.3, 5)
    diffs = [fc[i + 1] - fc[i] for i in range(len(fc) - 1)]
    assert max(diffs) - min(diffs) < 1e-9


def test_holt_winters_recovers_seasonal_pattern():
    pat = [10, -5, 3, -8]
    y = [pat[i % 4] for i in range(40)]
    _, _, _, fc = holt_winters_add(y, 0.3, 0.05, 0.3, 4, 8)
    assert max(abs(fc[i] - pat[i % 4]) for i in range(8)) < 2.0


def test_holt_winters_recovers_trend_with_season():
    pat = [10, -5, 3, -8]
    y = [pat[i % 4] + 0.5 * i for i in range(40)]
    _, trend, _, _ = holt_winters_add(y, 0.3, 0.1, 0.3, 4, 4)
    assert abs(trend - 0.5) < 0.15


def test_validation():
    with pytest.raises(ValueError):
        holt_linear([1.0], 0.5, 0.5)             # too short
    with pytest.raises(ValueError):
        holt_linear([1.0, 2.0], 1.5, 0.5)        # alpha out of range
    with pytest.raises(ValueError):
        holt_winters_add([1, 2, 3], 0.5, 0.5, 0.5, 4)   # < 2 seasons
    with pytest.raises(ValueError):
        holt_winters_add([float(i) for i in range(40)], 0.5, 0.5, 0.5, 1)  # period < 2
