"""Kruskal-Wallis H and Friedman rank tests."""

import math
import random

import pytest

from quantforge import kruskal_wallis_test, friedman_test, mann_whitney_u


def test_kruskal_wallis_textbook():
    g1 = [2.9, 3.0, 2.5, 2.6, 3.2]
    g2 = [3.8, 2.7, 4.0, 2.4]
    g3 = [2.8, 3.4, 3.7, 2.2, 2.0]
    r = kruskal_wallis_test(g1, g2, g3)
    assert abs(r["statistic"] - 0.7714) < 1e-3
    assert r["df"] == 2
    assert abs(r["p_value"] - 0.68) < 0.01


def test_two_group_equals_mann_whitney_z_squared():
    a = random.Random(1).sample(range(1000), 25)
    b = random.Random(2).sample(range(1000, 2000), 25)
    h = kruskal_wallis_test(a, b)["statistic"]
    u, _ = mann_whitney_u(a, b)
    n1 = n2 = 25
    N = 50
    mean = n1 * n2 / 2
    var = n1 * n2 * (N + 1) / 12
    z = (u - mean) / math.sqrt(var)
    assert abs(h - z * z) < 1e-6


def test_kruskal_wallis_detects_difference():
    rng = random.Random(7)
    gs = [[rng.gauss(0, 1) for _ in range(30)],
          [rng.gauss(1, 1) for _ in range(30)],
          [rng.gauss(2, 1) for _ in range(30)]]
    assert kruskal_wallis_test(*gs)["p_value"] < 0.001


def test_kruskal_wallis_null_calibration():
    rng = random.Random(5)
    rej = 0
    N = 1000
    for _ in range(N):
        gs = [[rng.gauss(0, 1) for _ in range(15)] for _ in range(3)]
        if kruskal_wallis_test(*gs)["p_value"] < 0.05:
            rej += 1
    assert 0.03 <= rej / N <= 0.07


def test_friedman_textbook():
    data = [[1, 2, 3], [2, 3, 1], [3, 1, 2], [1, 2, 3], [2, 3, 1]]
    r = friedman_test(data)
    assert abs(r["statistic"] - 0.4) < 1e-9
    assert r["df"] == 2


def test_friedman_detects_consistent_ordering():
    rng = random.Random(9)
    blocks = []
    for _ in range(12):
        base = rng.gauss(0, 1)
        blocks.append([base, base + 1.0, base + 2.0])   # treatment 3 always best
    r = friedman_test(blocks)
    assert r["p_value"] < 0.001


def test_friedman_null_calibration():
    rng = random.Random(11)
    rej = 0
    N = 1000
    for _ in range(N):
        blocks = [[rng.gauss(0, 1) for _ in range(4)] for _ in range(10)]
        if friedman_test(blocks)["p_value"] < 0.05:
            rej += 1
    assert 0.03 <= rej / N <= 0.07


def test_validation():
    with pytest.raises(ValueError):
        kruskal_wallis_test([1, 2, 3])                 # only one group
    with pytest.raises(ValueError):
        friedman_test([[1, 2, 3]])                     # only one block
    with pytest.raises(ValueError):
        friedman_test([[1, 2], [1, 2, 3]])             # ragged
