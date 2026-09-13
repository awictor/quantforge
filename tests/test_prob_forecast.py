"""Probabilistic-forecast scoring rules."""

import random
import statistics

import pytest

from quantforge import pinball_loss, interval_score, interval_coverage, crps_ensemble


def test_pinball_minimized_at_true_quantile():
    rng = random.Random(1)
    data = [rng.gauss(0, 1) for _ in range(20000)]
    tau = 0.9
    truth = statistics.quantiles(data, n=100)[88]
    at_truth = pinball_loss(data, [truth] * len(data), tau)
    below = pinball_loss(data, [truth - 0.5] * len(data), tau)
    above = pinball_loss(data, [truth + 0.5] * len(data), tau)
    assert at_truth < below
    assert at_truth < above


def test_pinball_asymmetry():
    # tau=0.9 penalizes under-forecasts (a > q) more than over-forecasts.
    over = pinball_loss([10.0], [12.0], 0.9)     # forecast too high by 2
    under = pinball_loss([10.0], [8.0], 0.9)     # forecast too low by 2
    assert under > over


def test_crps_single_member_is_abs_error():
    assert crps_ensemble(3.0, [5.0]) == 2.0


def test_crps_perfect_is_zero():
    assert crps_ensemble(4.0, [4.0, 4.0, 4.0]) == 0.0


def test_crps_non_negative():
    rng = random.Random(2)
    ens = [rng.gauss(0, 1) for _ in range(50)]
    assert crps_ensemble(0.3, ens) >= 0.0


def test_coverage_counts_inside():
    a = [0, 1, 2, 3, 4]
    lo = [-1, -1, -1, -1, -1]
    up = [1, 1, 3, 3, 3]
    assert interval_coverage(a, lo, up) == 0.8


def test_interval_score_rewards_tight_covering():
    rng = random.Random(3)
    act = [rng.gauss(0, 1) for _ in range(1000)]
    tight = interval_score(act, [-1.65] * 1000, [1.65] * 1000, alpha=0.1)
    wide = interval_score(act, [-5.0] * 1000, [5.0] * 1000, alpha=0.1)
    assert tight < wide


def test_interval_score_penalizes_misses():
    # An interval that never covers scores worse than one that always covers.
    act = [0.0] * 100
    covering = interval_score(act, [-1.0] * 100, [1.0] * 100, alpha=0.1)
    missing = interval_score(act, [2.0] * 100, [3.0] * 100, alpha=0.1)
    assert missing > covering


def test_validation():
    with pytest.raises(ValueError):
        pinball_loss([1.0], [1.0], tau=1.0)
    with pytest.raises(ValueError):
        interval_score([1.0], [0.0], [2.0], alpha=0.0)
    with pytest.raises(ValueError):
        crps_ensemble(1.0, [])
    with pytest.raises(ValueError):
        interval_coverage([1.0], [0.0], [])
