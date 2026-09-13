"""Wald-Wolfowitz runs test."""

import random

import pytest

from quantforge import runs_test, runs_test_binary


def test_random_holds_size():
    rej = 0
    trials = 2000
    for s in range(trials):
        r = random.Random(s)
        x = [r.gauss(0, 1) for _ in range(200)]
        _, p = runs_test(x)
        if p < 0.05:
            rej += 1
    assert 0.03 < rej / trials < 0.08


def test_sorted_data_too_few_runs():
    rng = random.Random(1)
    x = sorted(rng.gauss(0, 1) for _ in range(200))
    z, p = runs_test(x)
    assert z < 0
    assert p < 0.001


def test_alternating_too_many_runs():
    z, p = runs_test([1, -1] * 100)
    assert z > 0
    assert p < 0.001


def test_binary_known_value():
    # 1,1,1,0,0,0 -> 2 runs, n1=n2=3; expected runs = 4, so z is negative.
    z, _ = runs_test_binary([1, 1, 1, 0, 0, 0])
    assert z < 0
    assert abs(z + 1.8257) < 1e-3


def test_binary_alternating():
    z, p = runs_test_binary([0, 1] * 50)
    assert z > 0
    assert p < 0.001


def test_validation():
    with pytest.raises(ValueError):
        runs_test([1.0])
    with pytest.raises(ValueError):
        runs_test_binary([1, 1, 1])          # only one distinct symbol
    with pytest.raises(ValueError):
        runs_test_binary([1, 2, 3])          # three symbols
