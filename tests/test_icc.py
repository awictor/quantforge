"""Intraclass correlation coefficient (Shrout-Fleiss)."""

import random

import pytest

from quantforge import icc


def test_shrout_fleiss_canonical():
    # The canonical table from Shrout & Fleiss (1979).
    data = [[9, 2, 5, 8], [6, 1, 3, 2], [8, 4, 6, 8],
            [7, 1, 2, 6], [10, 5, 6, 9], [6, 2, 4, 7]]
    r = icc(data)
    assert abs(r["icc1"] - 0.166) < 0.005
    assert abs(r["icc2_1"] - 0.290) < 0.005
    assert abs(r["icc3_1"] - 0.715) < 0.005


def test_perfect_agreement():
    r = icc([[5, 5], [3, 3], [8, 8], [1, 1]])
    for v in r.values():
        assert abs(v - 1.0) < 1e-9


def test_noise_near_zero():
    rng = random.Random(3)
    data = [[rng.gauss(0, 1) for _ in range(4)] for _ in range(100)]
    r = icc(data)
    assert abs(r["icc2_1"]) < 0.1
    assert abs(r["icc3_1"]) < 0.1


def test_spearman_brown_relation():
    data = [[9, 2, 5, 8], [6, 1, 3, 2], [8, 4, 6, 8],
            [7, 1, 2, 6], [10, 5, 6, 9], [6, 2, 4, 7]]
    r = icc(data)
    k = 4
    sb = k * r["icc3_1"] / (1 + (k - 1) * r["icc3_1"])
    assert abs(r["icc3_k"] - sb) < 1e-9
    sb1 = k * r["icc1"] / (1 + (k - 1) * r["icc1"])
    assert abs(r["icc1k"] - sb1) < 1e-9


def test_consistency_exceeds_agreement_under_offset():
    # rater 2 = rater 1 + 2: perfect consistency, imperfect absolute agreement.
    data = [[10, 12], [20, 22], [30, 32], [40, 42], [50, 52]]
    r = icc(data)
    assert abs(r["icc3_1"] - 1.0) < 1e-9
    assert r["icc3_1"] > r["icc2_1"]


def test_validation():
    with pytest.raises(ValueError):
        icc([[1, 2]])                      # one subject
    with pytest.raises(ValueError):
        icc([[1], [2]])                    # one rater
    with pytest.raises(ValueError):
        icc([[1, 2], [3, 4, 5]])           # ragged
