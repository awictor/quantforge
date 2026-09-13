"""Nonparametric two-sample scale tests: Ansari-Bradley and Mood."""

import random

import pytest

from quantforge import ansari_bradley_test, mood_test


def test_ansari_tiny_statistic():
    # x=[1,2] ranks 1,2 in pool of 4 -> scores min(r,5-r)=1,2 -> sum 3.
    assert ansari_bradley_test([1, 2], [3, 4])["statistic"] == 3.0


def test_ansari_null_calibration():
    rng = random.Random(3)
    zs = [ansari_bradley_test([rng.gauss(0, 1) for _ in range(30)],
                              [rng.gauss(0, 1) for _ in range(30)])["z"]
          for _ in range(500)]
    mean = sum(zs) / len(zs)
    var = sum(z * z for z in zs) / len(zs)
    assert abs(mean) < 0.2
    assert 0.8 < var < 1.25


def test_ansari_detects_scale():
    rng = random.Random(7)
    x = [rng.gauss(0, 1) for _ in range(80)]
    y = [rng.gauss(0, 4) for _ in range(80)]
    assert ansari_bradley_test(x, y)["p_value"] < 0.01


def test_mood_null_calibration():
    rng = random.Random(5)
    zs = [mood_test([rng.gauss(0, 1) for _ in range(30)],
                    [rng.gauss(0, 1) for _ in range(30)])["z"]
          for _ in range(500)]
    mean = sum(zs) / len(zs)
    var = sum(z * z for z in zs) / len(zs)
    assert abs(mean) < 0.2
    assert 0.8 < var < 1.25


def test_mood_detects_scale():
    rng = random.Random(9)
    x = [rng.gauss(0, 1) for _ in range(80)]
    y = [rng.gauss(0, 4) for _ in range(80)]
    assert mood_test(x, y)["p_value"] < 0.01


def test_rank_sum_blind_to_scale_but_these_are_not():
    # Same median, very different spread: a location test would not flag it.
    rng = random.Random(21)
    x = [rng.gauss(5, 1) for _ in range(100)]
    y = [rng.gauss(5, 5) for _ in range(100)]
    assert ansari_bradley_test(x, y)["p_value"] < 0.05
    assert mood_test(x, y)["p_value"] < 0.05


def test_validation():
    with pytest.raises(ValueError):
        ansari_bradley_test([], [1, 2])
    with pytest.raises(ValueError):
        mood_test([1, 2], [])
