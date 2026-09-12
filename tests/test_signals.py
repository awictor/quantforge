"""Trend and momentum signals."""

import pytest

from quantforge import sma, ema, macd, rsi, rolling_zscore, time_series_momentum


def test_sma():
    assert sma([1, 2, 3, 4, 5], 3) == [2, 3, 4]


def test_ema_reacts_faster_than_sma_after_step():
    series = [10] * 10 + [20] * 5
    e = ema(series, 5)
    s = sma(series, 5)
    # First post-step points: EMA above the corresponding lagging SMA window.
    assert e[10] > s[6]


def test_macd_positive_on_uptrend_and_definition():
    up = [float(i) for i in range(50)]
    ml, sl, hist = macd(up)
    assert ml[-1] > 0
    assert ml[-1] == pytest.approx(ema(up, 12)[-1] - ema(up, 26)[-1])


def test_rsi_range_and_trend():
    up = [float(i) for i in range(50)]
    down = [float(50 - i) for i in range(50)]
    ru, rd = rsi(up, 14), rsi(down, 14)
    assert all(0 <= x <= 100 for x in ru)
    assert ru[-1] > 70
    assert rd[-1] < 30


def test_rolling_zscore_flags_spike():
    z = rolling_zscore([1, 1, 1, 1, 10], 5)
    assert z[0] > 1.5


def test_time_series_momentum_sign():
    assert time_series_momentum([1, 2, 3, 2, 1], 2) == [1, 0, -1]


def test_validation():
    with pytest.raises(ValueError):
        sma([1, 2], 5)
    with pytest.raises(ValueError):
        macd([1.0] * 30, fast=26, slow=12)   # fast >= slow
