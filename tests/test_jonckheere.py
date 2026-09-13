"""Jonckheere-Terpstra ordered-trend test."""

import random

import pytest

from quantforge import jonckheere_terpstra_test


def test_null_mean_and_variance_formula():
    # 3 groups of 4, no ties: mean=(144-48)/4=24, var=(12*11*29 - 3*4*3*13)/72.
    r = jonckheere_terpstra_test([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
    assert abs(r["mean"] - 24.0) < 1e-9
    assert abs(r["variance"] - 3360.0 / 72.0) < 1e-9


def test_perfect_increasing_max_statistic():
    r = jonckheere_terpstra_test([[1, 2], [3, 4], [5, 6]])
    assert r["statistic"] == 12.0        # 3 pairs x 4 concordant each
    assert r["z"] > 0


def test_increasing_vs_decreasing_sign():
    inc = jonckheere_terpstra_test([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
    dec = jonckheere_terpstra_test([[9, 10, 11, 12], [5, 6, 7, 8], [1, 2, 3, 4]])
    assert inc["z"] > 3
    assert abs(inc["z"] + dec["z"]) < 1e-9        # mirror images


def test_null_calibration():
    rng = random.Random(5)
    zs = [jonckheere_terpstra_test([[rng.gauss(0, 1) for _ in range(8)]
                                    for _ in range(4)])["z"]
          for _ in range(2000)]
    mean = sum(zs) / len(zs)
    var = sum(z * z for z in zs) / len(zs)
    assert abs(mean) < 0.15
    assert 0.85 < var < 1.2


def test_detects_monotone_trend():
    rng = random.Random(7)
    gs = [[rng.gauss(m, 1) for _ in range(20)] for m in (0.0, 0.5, 1.0, 1.5)]
    assert jonckheere_terpstra_test(gs)["p_value"] < 0.01


def test_ties_handled():
    gs = [[1, 1, 2], [2, 2, 3], [3, 3, 4]]
    r = jonckheere_terpstra_test(gs)
    assert r["variance"] > 0
    assert 0.0 <= r["p_value"] <= 1.0


def test_validation():
    with pytest.raises(ValueError):
        jonckheere_terpstra_test([[1, 2, 3]])
