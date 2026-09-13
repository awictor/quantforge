"""Binomial-proportion confidence intervals."""

import random

import pytest

from quantforge import (wald_interval, wilson_interval, agresti_coull_interval,
                        clopper_pearson_interval)

ALL = [wald_interval, wilson_interval, agresti_coull_interval,
       clopper_pearson_interval]


def test_clopper_pearson_textbook():
    lo, hi = clopper_pearson_interval(2, 10, 0.95)
    assert abs(lo - 0.02521073) < 1e-5
    assert abs(hi - 0.55609546) < 1e-5


def test_all_contain_point_estimate():
    for f in ALL:
        lo, hi = f(2, 10)
        assert lo <= 0.2 <= hi
        assert 0.0 <= lo <= hi <= 1.0


def test_wilson_stays_in_unit_interval_at_extremes():
    lo, hi = wilson_interval(0, 20)
    assert lo == 0.0 and 0.0 < hi < 1.0
    lo, hi = wilson_interval(20, 20)
    assert hi == 1.0 and 0.0 < lo < 1.0


def test_clopper_pearson_one_sided_at_edges():
    lo, hi = clopper_pearson_interval(0, 10)
    assert lo == 0.0 and hi > 0.0
    lo, hi = clopper_pearson_interval(10, 10)
    assert hi == 1.0 and lo < 1.0


def test_wald_approaches_wilson_for_large_n():
    lo1, hi1 = wald_interval(500, 1000)
    lo2, hi2 = wilson_interval(500, 1000)
    assert abs(lo1 - lo2) < 1e-2
    assert abs(hi1 - hi2) < 1e-2


def test_higher_confidence_widens():
    for f in ALL:
        w90 = f(8, 25, 0.90)
        w99 = f(8, 25, 0.99)
        assert (w99[1] - w99[0]) > (w90[1] - w90[0])


def test_clopper_pearson_contains_wilson_width():
    # The exact interval is the widest (most conservative) of the four here.
    k, n = 3, 20
    cp = clopper_pearson_interval(k, n)
    wil = wilson_interval(k, n)
    assert (cp[1] - cp[0]) >= (wil[1] - wil[0]) - 1e-9


def test_clopper_pearson_coverage_at_least_nominal():
    p, n = 0.3, 40
    rng = random.Random(11)
    hits, trials = 0, 4000
    for _ in range(trials):
        k = sum(1 for _ in range(n) if rng.random() < p)
        lo, hi = clopper_pearson_interval(k, n)
        if lo <= p <= hi:
            hits += 1
    assert hits / trials >= 0.95


def test_validation():
    for f in ALL:
        with pytest.raises(ValueError):
            f(11, 10)
        with pytest.raises(ValueError):
            f(0, 0)
        with pytest.raises(ValueError):
            f(2, 10, confidence=1.0)
