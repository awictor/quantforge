"""Bollinger bands, ATR, and Donchian channel."""

import random

import pytest

from quantforge import (
    bollinger_bands, average_true_range, donchian_channel, sma,
)


def _prices():
    random.seed(4)
    return [100 + random.gauss(0, 2) for _ in range(100)]


def test_bollinger_ordered_and_middle_is_sma():
    prices = _prices()
    lo, mid, up = bollinger_bands(prices, 20, 2.0)
    assert all(lo[i] <= mid[i] <= up[i] for i in range(len(mid)))
    assert all(mid[i] == pytest.approx(sma(prices, 20)[i]) for i in range(len(mid)))


def test_bollinger_coverage():
    prices = _prices()
    lo, mid, up = bollinger_bands(prices, 20, 2.0)
    inside = sum(1 for i in range(len(mid)) if lo[i] <= prices[i + 19] <= up[i])
    assert inside / len(mid) > 0.85


def test_atr_positive():
    prices = _prices()
    highs = [p + 1 for p in prices]
    lows = [p - 1 for p in prices]
    assert all(x > 0 for x in average_true_range(highs, lows, prices, 14))


def test_donchian_ordered_and_max():
    prices = _prices()
    highs = [p + 1 for p in prices]
    lows = [p - 1 for p in prices]
    lo, up = donchian_channel(highs, lows, 20)
    assert all(up[i] >= lo[i] for i in range(len(up)))
    assert up[-1] == pytest.approx(max(highs[-20:]))


def test_validation():
    prices = _prices()
    with pytest.raises(ValueError):
        bollinger_bands(prices, 1)
    with pytest.raises(ValueError):
        average_true_range(prices[:5], prices, prices, 14)
