"""Dunn's post-hoc pairwise test."""

import random

import pytest

from quantforge import dunn_test, kruskal_wallis_test


def test_pair_count():
    res = dunn_test([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert len(res) == 3
    assert {r["groups"] for r in res} == {(0, 1), (0, 2), (1, 2)}


def test_separated_groups_significant():
    res = dunn_test([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15]])
    far = next(r for r in res if r["groups"] == (0, 2))
    assert far["p_adjusted"] < 0.05


def test_identical_groups_null():
    res = dunn_test([[1, 2, 3, 4, 5], [1, 2, 3, 4, 5]], adjust=None)
    assert abs(res[0]["z"]) < 1e-9
    assert abs(res[0]["p_value"] - 1.0) < 1e-9


def test_z_sign():
    res = dunn_test([[1, 2, 3], [10, 11, 12]], adjust=None)
    assert res[0]["z"] < 0            # group 0 has the lower ranks


def test_adjusted_at_least_raw():
    res = dunn_test([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15]])
    assert all(r["p_adjusted"] >= r["p_value"] - 1e-12 for r in res)


def test_z_squared_equals_kruskal_h_for_two_groups():
    a = random.Random(1).sample(range(1000), 20)
    b = random.Random(2).sample(range(1000, 2000), 20)
    res = dunn_test([a, b], adjust=None)
    h = kruskal_wallis_test(a, b)["statistic"]
    assert abs(res[0]["z"] ** 2 - h) < 1e-6


def test_bonferroni_adjust():
    res = dunn_test([[1, 2, 3], [4, 5, 6], [7, 8, 9]], adjust="bonferroni")
    for r in res:
        assert r["p_adjusted"] <= min(1.0, r["p_value"] * 3) + 1e-12


def test_validation():
    with pytest.raises(ValueError):
        dunn_test([[1, 2, 3]])
    with pytest.raises(ValueError):
        dunn_test([[1, 2], [3, 4]], adjust="bad")
