"""Anderson-Darling normality test."""

import random

import pytest

from quantforge import anderson_darling_normal


def test_normal_holds_size():
    rej = 0
    trials = 2000
    for s in range(trials):
        r = random.Random(s)
        x = [r.gauss(0, 1) for _ in range(100)]
        _, p = anderson_darling_normal(x)
        if p < 0.05:
            rej += 1
    assert 0.03 < rej / trials < 0.08


def test_normal_not_rejected():
    r = random.Random(1)
    x = [r.gauss(5, 2) for _ in range(500)]
    _, p = anderson_darling_normal(x)
    assert p > 0.05


def test_exponential_rejected():
    r = random.Random(2)
    x = [r.expovariate(1.0) for _ in range(500)]
    a, p = anderson_darling_normal(x)
    assert a > 1.0
    assert p < 0.001


def test_heavy_tails_rejected():
    r = random.Random(3)
    x = [r.gauss(0, 1) if r.random() < 0.9 else r.gauss(0, 8) for _ in range(1000)]
    _, p = anderson_darling_normal(x)
    assert p < 0.001


def test_statistic_non_negative():
    r = random.Random(4)
    x = [r.gauss(0, 1) for _ in range(50)]
    a, _ = anderson_darling_normal(x)
    assert a >= 0.0


def test_validation():
    with pytest.raises(ValueError):
        anderson_darling_normal([1.0, 2.0, 3.0])          # < 8
    with pytest.raises(ValueError):
        anderson_darling_normal([2.0] * 10)               # zero variance
