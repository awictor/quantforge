"""Median and restricted-mean survival from Kaplan-Meier."""

import math
import random

import pytest

from quantforge import median_survival_time, restricted_mean_survival_time


def test_median_no_censoring():
    t = list(range(1, 11))
    e = [1] * 10
    assert median_survival_time(t, e) == 5


def test_median_none_when_all_censored():
    assert median_survival_time([1, 2, 3], [0, 0, 0]) is None


def test_rmst_by_hand():
    # times 1,2,3 all events: S = 1, 2/3, 1/3 on the unit steps; RMST(3) = 2.
    assert abs(restricted_mean_survival_time([1, 2, 3], [1, 1, 1], 3) - 2.0) < 1e-9


def test_rmst_exponential():
    lam, tau = 0.5, 3.0
    rng = random.Random(1)
    t = [rng.expovariate(lam) for _ in range(20000)]
    e = [1] * 20000
    approx = restricted_mean_survival_time(t, e, tau)
    exact = (1 - math.exp(-lam * tau)) / lam
    assert abs(approx - exact) < 0.05


def test_rmst_increases_with_horizon():
    t = list(range(1, 21))
    e = [1] * 20
    assert restricted_mean_survival_time(t, e, 5) < restricted_mean_survival_time(t, e, 15)


def test_rmst_capped_by_tau():
    # RMST can never exceed the horizon.
    t = [100.0] * 10
    e = [1] * 10
    assert restricted_mean_survival_time(t, e, 5.0) <= 5.0 + 1e-9


def test_validation():
    with pytest.raises(ValueError):
        restricted_mean_survival_time([1, 2], [1, 1], tau=0.0)
    with pytest.raises(ValueError):
        median_survival_time([], [])
