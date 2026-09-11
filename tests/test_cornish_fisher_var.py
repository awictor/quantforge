"""Cornish-Fisher (skew/kurtosis-adjusted) VaR (perfmetrics module)."""

import math
import random

import pytest

from quantforge import cornish_fisher_var


def test_reduces_to_parametric_for_normal():
    rng = random.Random(3)
    norm = [rng.gauss(0.0, 0.01) for _ in range(20000)]
    mu = sum(norm) / len(norm)
    sd = (sum((x - mu) ** 2 for x in norm) / (len(norm) - 1)) ** 0.5
    z = 1.6448536269514722
    param = -(mu - z * sd)
    assert cornish_fisher_var(norm) == pytest.approx(param, abs=1e-3)


def test_positive_loss():
    rng = random.Random(5)
    x = [rng.gauss(0.0, 0.01) for _ in range(5000)]
    assert cornish_fisher_var(x) > 0.0


def test_negative_skew_increases_var():
    rng = random.Random(7)
    base = [rng.gauss(0.0, 0.01) for _ in range(5000)]
    skewed = base + [-0.1, -0.12, -0.15]   # fat left tail
    mu = sum(skewed) / len(skewed)
    sd = (sum((x - mu) ** 2 for x in skewed) / (len(skewed) - 1)) ** 0.5
    z = 1.6448536269514722
    gaussian = -(mu - z * sd)
    assert cornish_fisher_var(skewed) > gaussian


def test_higher_confidence_larger():
    rng = random.Random(9)
    x = [rng.gauss(0.0, 0.01) for _ in range(5000)]
    assert cornish_fisher_var(x, 0.99) > cornish_fisher_var(x, 0.95)


def test_horizon_scaling_between_linear_and_sqrt():
    rng = random.Random(11)
    x = [rng.gauss(0.001, 0.01) for _ in range(5000)]
    base = cornish_fisher_var(x, 0.95, 1.0)
    quad = cornish_fisher_var(x, 0.95, 4.0)
    # Deviation scales with sqrt(4)=2, mean linearly with 4: VaR between.
    assert base < quad < 5 * base


def test_validation():
    with pytest.raises(ValueError):
        cornish_fisher_var([0.01, 0.01])   # zero variance
    with pytest.raises(ValueError):
        cornish_fisher_var([0.01])         # need >= 2
