"""Engle ARCH-LM test for conditional heteroskedasticity."""

import math
import random

import pytest

from quantforge import arch_lm_test


def test_iid_holds_size():
    rej = 0
    trials = 2000
    for s in range(trials):
        r = random.Random(s)
        x = [r.gauss(0, 1) for _ in range(300)]
        _, p = arch_lm_test(x, lags=2)
        if p < 0.05:
            rej += 1
    assert 0.03 < rej / trials < 0.08


def test_garch_rejected():
    r = random.Random(1)
    omega, alpha, beta = 0.05, 0.15, 0.80
    h = omega / (1 - alpha - beta)
    x = []
    for _ in range(2000):
        z = r.gauss(0, 1)
        e = math.sqrt(h) * z
        x.append(e)
        h = omega + alpha * e * e + beta * h
    lm, p = arch_lm_test(x, lags=5)
    assert lm > 20.0
    assert p < 1e-4


def test_single_iid_not_rejected():
    r = random.Random(3)
    y = [r.gauss(0, 2) for _ in range(500)]
    _, p = arch_lm_test(y, lags=3)
    assert p > 0.01


def test_statistic_non_negative():
    r = random.Random(4)
    x = [r.gauss(0, 1) for _ in range(100)]
    lm, _ = arch_lm_test(x, lags=1)
    assert lm >= 0.0


def test_validation():
    with pytest.raises(ValueError):
        arch_lm_test([1.0, 2.0, 3.0], lags=0)
    with pytest.raises(ValueError):
        arch_lm_test([1.0, 2.0, 3.0], lags=5)      # too short
