"""Log-rank (Mantel-Cox) two-sample survival test."""

import random

import pytest

from quantforge import log_rank_test


def _sample(rate, n, seed):
    r = random.Random(seed)
    t, e = [], []
    for _ in range(n):
        x = r.expovariate(rate)
        c = r.expovariate(rate * 0.3)     # independent censoring
        if x <= c:
            t.append(x)
            e.append(1)
        else:
            t.append(c)
            e.append(0)
    return t, e


def test_identical_groups_not_rejected():
    t1, e1 = _sample(1.0, 200, 1)
    t2, e2 = _sample(1.0, 200, 2)
    _, p = log_rank_test(t1, e1, t2, e2)
    assert p > 0.05


def test_different_hazards_rejected():
    t1, e1 = _sample(1.0, 200, 3)
    t2, e2 = _sample(3.0, 200, 4)      # much faster hazard
    chi2, p = log_rank_test(t1, e1, t2, e2)
    assert chi2 > 10.0
    assert p < 0.001


def test_same_data_zero_statistic():
    t, e = _sample(1.0, 100, 5)
    chi2, p = log_rank_test(t, e, t, e)
    assert abs(chi2) < 1e-9
    assert abs(p - 1.0) < 1e-9


def test_symmetry():
    t1, e1 = _sample(1.0, 150, 6)
    t2, e2 = _sample(2.0, 150, 7)
    a, _ = log_rank_test(t1, e1, t2, e2)
    b, _ = log_rank_test(t2, e2, t1, e1)
    assert abs(a - b) < 1e-9        # (O-E)^2 is direction-invariant


def test_validation():
    with pytest.raises(ValueError):
        log_rank_test([], [], [1.0], [1])
