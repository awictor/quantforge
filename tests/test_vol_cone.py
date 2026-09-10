"""Tests for the realized-volatility cone."""

import math
import random

import pytest

from quantforge import vol_cone, close_to_close


def _gbm_closes(n=500, sigma=0.2, seed=1):
    rng = random.Random(seed)
    closes = [100.0]
    for _ in range(n):
        closes.append(closes[-1] * math.exp(-0.5 * sigma * sigma / 252
                                            + sigma / math.sqrt(252) * rng.gauss(0, 1)))
    return closes


def test_percentiles_are_ordered():
    cone = vol_cone(_gbm_closes(), [5, 21, 63, 126])
    for p in cone:
        assert p.minimum <= p.p25 <= p.median <= p.p75 <= p.maximum


def test_current_within_min_max():
    cone = vol_cone(_gbm_closes(), [21, 63])
    for p in cone:
        assert p.minimum <= p.current <= p.maximum


def test_median_near_true_vol():
    cone = vol_cone(_gbm_closes(sigma=0.2), [21, 63, 126])
    for p in cone:
        assert p.median == pytest.approx(0.2, abs=0.05)


def test_cone_narrows_with_longer_windows():
    # The spread of realized vol shrinks as the estimation window grows.
    cone = vol_cone(_gbm_closes(), [5, 21, 63, 126])
    spreads = [p.maximum - p.minimum for p in cone]
    assert spreads == sorted(spreads, reverse=True)


def test_one_window_equals_close_to_close():
    # A single window spanning the whole return series matches close_to_close.
    closes = _gbm_closes(n=100)          # 100 appends + seed = 101 closes -> 100 returns
    cone = vol_cone(closes, [100])       # one block of all 100 returns
    assert len(cone) == 1
    assert cone[0].current == pytest.approx(close_to_close(closes), abs=1e-9)


def test_skips_windows_too_long():
    closes = _gbm_closes(n=30)
    cone = vol_cone(closes, [5, 21, 1000])   # 1000 too long
    assert [p.window for p in cone] == [5, 21]


def test_windows_are_sorted_in_output():
    cone = vol_cone(_gbm_closes(), [63, 5, 21])
    assert [p.window for p in cone] == [5, 21, 63]
