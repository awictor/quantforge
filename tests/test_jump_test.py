"""Barndorff-Nielsen-Shephard realized-volatility jump test."""

import math
import random

import pytest

from quantforge import bns_jump_test, tripower_quarticity


def _path(sigma, n, dt, jump=0.0, seed=0):
    rng = random.Random(seed)
    r = [sigma * math.sqrt(dt) * rng.gauss(0, 1) for _ in range(n)]
    if jump:
        r[n // 2] += jump
    return r


def test_null_size_is_about_five_percent():
    sigma, n = 0.01, 1000
    rej = 0
    trials = 2000
    zs = []
    for s in range(trials):
        z, p = bns_jump_test(_path(sigma, n, 1.0 / n, 0.0, s))
        zs.append(z)
        if p < 0.05:
            rej += 1
    # Statistic is roughly N(0,1) and the test holds its size.
    assert abs(sum(zs) / trials) < 0.1
    assert 0.03 < rej / trials < 0.08


def test_high_power_against_a_jump():
    sigma, n, jump = 0.01, 1000, 0.06
    rej = 0
    trials = 400
    for s in range(trials):
        _, p = bns_jump_test(_path(sigma, n, 1.0 / n, jump, s))
        if p < 0.05:
            rej += 1
    assert rej / trials > 0.9


def test_jump_gives_large_positive_z():
    z, p = bns_jump_test(_path(0.01, 1000, 1.0 / 1000, 0.06, 3))
    assert z > 3.0
    assert p < 0.01


def test_tripower_quarticity_consistent():
    sigma, n = 0.01, 3000
    tq = [tripower_quarticity(_path(sigma, n, 1.0 / n, 0.0, s)) for s in range(40)]
    mean = sum(tq) / 40
    assert abs(mean - sigma ** 4) < 0.2 * sigma ** 4


def test_tripower_quarticity_jump_robust():
    # A single jump barely moves the tripower quarticity relative to RV^-based ones.
    clean = tripower_quarticity(_path(0.01, 2000, 1.0 / 2000, 0.0, 5))
    jumped = tripower_quarticity(_path(0.01, 2000, 1.0 / 2000, 0.06, 5))
    assert jumped < 5.0 * clean


def test_validation():
    with pytest.raises(ValueError):
        bns_jump_test([0.01, 0.02])
    with pytest.raises(ValueError):
        tripower_quarticity([0.01, 0.02])
