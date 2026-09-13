"""Streaming exponentially-weighted mean and variance."""

import random

import pytest

from quantforge import EWMAStats, ewma


def test_constant_series():
    e = EWMAStats(0.9, [5.0] * 100)
    assert e.mean == 5.0
    assert e.variance() == 0.0


def test_batch_matches_manual_recursion():
    lam = 0.94
    vals = [1, 2, 3, 4, 5]
    manual = [float(vals[0])]
    for i in range(1, 5):
        manual.append(lam * manual[-1] + (1 - lam) * vals[i])
    assert ewma(vals, lam) == manual


def test_tracks_step_change():
    e = EWMAStats(0.9)
    for _ in range(50):
        e.update(0.0)
    for _ in range(50):
        e.update(10.0)
    assert 9.0 < e.mean < 10.0        # moved toward the new level, lagging slightly


def test_vol_on_normal():
    rng = random.Random(3)
    e = EWMAStats(0.97)
    for _ in range(50000):
        e.update(rng.gauss(0, 2))
    assert 1.5 < e.std() < 2.3


def test_reacts_to_vol_regime_change():
    rng = random.Random(5)
    e = EWMAStats(0.9)
    for _ in range(200):
        e.update(rng.gauss(0, 0.5))
    low = e.std()
    for _ in range(200):
        e.update(rng.gauss(0, 5))
    high = e.std()
    assert high > low * 3


def test_first_point_seeds_mean():
    e = EWMAStats(0.9)
    e.update(7.0)
    assert e.mean == 7.0
    assert e.variance() == 0.0


def test_validation():
    with pytest.raises(ValueError):
        EWMAStats(1.5)
    with pytest.raises(ValueError):
        ewma([])
