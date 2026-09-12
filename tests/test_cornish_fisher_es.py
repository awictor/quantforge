"""Cornish-Fisher (skew/kurtosis-adjusted) expected shortfall."""

import math
import random

import pytest

from quantforge import (cornish_fisher_expected_shortfall, cornish_fisher_var)
from quantforge.mathfns import norm_ppf


def test_normal_reduces_to_gaussian_es():
    rng = random.Random(42)
    x = [rng.gauss(0.0, 0.02) for _ in range(200000)]
    c = 0.95
    mu = sum(x) / len(x)
    sd = (sum((v - mu) ** 2 for v in x) / (len(x) - 1)) ** 0.5
    z = norm_ppf(1 - c)
    phi = math.exp(-0.5 * z * z) / math.sqrt(2 * math.pi)
    gauss_es = -(mu - phi / (1 - c) * sd)
    assert abs(cornish_fisher_expected_shortfall(x, c) - gauss_es) < 5e-4


def test_es_at_least_var():
    rng = random.Random(1)
    x = [rng.gauss(0.001, 0.02) for _ in range(50000)]
    for c in (0.9, 0.95, 0.99):
        assert (cornish_fisher_expected_shortfall(x, c)
                >= cornish_fisher_var(x, c) - 1e-9)


def test_matches_empirical_under_mild_nonnormality():
    rng = random.Random(3)
    x = [0.001 + 0.02 * (z - 0.08 * (z * z - 1)) for z in
         (rng.gauss(0, 1) for _ in range(300000))]
    c = 0.95
    cf = cornish_fisher_expected_shortfall(x, c)
    s = sorted(x)
    k = int((1 - c) * len(s))
    emp = -sum(s[:k]) / k
    assert abs(cf - emp) / emp < 0.02


def test_horizon_scales_deviation():
    rng = random.Random(5)
    x = [rng.gauss(0.0, 0.02) for _ in range(20000)]
    one = cornish_fisher_expected_shortfall(x, 0.95, horizon=1.0)
    four = cornish_fisher_expected_shortfall(x, 0.95, horizon=4.0)
    # With ~zero mean the loss scales with sqrt(horizon).
    assert abs(four / one - 2.0) < 0.05


def test_fat_left_tail_raises_es():
    rng = random.Random(7)
    normal = [rng.gauss(0.0, 0.02) for _ in range(100000)]
    fat = []
    rng2 = random.Random(8)
    for _ in range(100000):
        if rng2.random() < 0.03:
            fat.append(rng2.gauss(-0.05, 0.03))
        else:
            fat.append(rng2.gauss(0.0015, 0.018))
    assert (cornish_fisher_expected_shortfall(fat, 0.95)
            > cornish_fisher_expected_shortfall(normal, 0.95))


def test_validation():
    with pytest.raises(ValueError):
        cornish_fisher_expected_shortfall([0.01], 0.95)
    with pytest.raises(ValueError):
        cornish_fisher_expected_shortfall([0.01, 0.01, 0.01], 0.95)  # zero variance
    with pytest.raises(ValueError):
        cornish_fisher_expected_shortfall([0.01, -0.02], 1.0)
