"""Generalized downside-risk ratios: Kappa and upside-potential."""

import math
import random

import pytest

from quantforge import kappa_ratio, upside_potential_ratio, lower_partial_moment
from quantforge.perfmetrics import sortino_ratio


def _returns(n, seed):
    rng = random.Random(seed)
    return [rng.gauss(0.001, 0.02) for _ in range(n)]


def test_kappa2_matches_sortino():
    r = _returns(500, 1)
    k2 = kappa_ratio(r, 0.0, 2)
    assert abs(k2 * math.sqrt(252) - sortino_ratio(r, 0.0, 0.0, 252)) < 1e-9


def test_kappa_decreases_with_order():
    r = _returns(500, 1)
    k1 = kappa_ratio(r, 0.0, 1)
    k2 = kappa_ratio(r, 0.0, 2)
    k3 = kappa_ratio(r, 0.0, 3)
    assert k1 > k2 > k3          # deeper shortfalls penalized more heavily


def test_higher_target_lowers_kappa():
    r = _returns(500, 1)
    assert kappa_ratio(r, 0.002, 2) < kappa_ratio(r, 0.0, 2)


def test_symmetric_kappa_near_zero():
    rng = random.Random(2)
    sym = [rng.gauss(0.0, 0.02) for _ in range(2000)]
    assert abs(kappa_ratio(sym, 0.0, 2)) < 0.1


def test_upside_potential_positive_for_positive_drift():
    r = _returns(500, 1)
    assert upside_potential_ratio(r, 0.0) > 0.0


def test_lower_partial_moment_nonnegative():
    r = _returns(200, 3)
    assert lower_partial_moment(r, 0.0, 2) >= 0.0


def test_validation():
    with pytest.raises(ValueError):
        kappa_ratio([0.01, 0.02, 0.03])          # no downside
    with pytest.raises(ValueError):
        kappa_ratio(_returns(50, 1), 0.0, 0)     # order <= 0
    with pytest.raises(ValueError):
        lower_partial_moment([], 0.0, 2)         # empty
