"""Probabilistic Sharpe ratio and minimum track record length."""

import random
import statistics

import pytest

from quantforge import (
    probabilistic_sharpe_ratio, minimum_track_record_length,
    deflated_sharpe_ratio,
)


def _sample():
    random.seed(5)
    return [random.gauss(0.001, 0.01) for _ in range(250)]


def test_psr_above_half_for_positive_edge():
    assert probabilistic_sharpe_ratio(_sample(), 0.0) > 0.5


def test_psr_lower_against_higher_benchmark():
    rets = _sample()
    assert probabilistic_sharpe_ratio(rets, 0.1) < probabilistic_sharpe_ratio(rets, 0.0)


def test_min_trl_positive():
    assert minimum_track_record_length(_sample(), 0.0, 0.95) > 0


def test_min_trl_longer_for_smaller_edge():
    rets = _sample()
    sr = statistics.mean(rets) / statistics.stdev(rets)
    assert minimum_track_record_length(rets, sr * 0.5, 0.95) > \
        minimum_track_record_length(rets, 0.0, 0.95)


def test_min_trl_requires_edge():
    with pytest.raises(ValueError):
        minimum_track_record_length(_sample(), 10.0, 0.95)


def test_dsr_single_trial_equals_psr():
    rets = _sample()
    assert deflated_sharpe_ratio(rets, 1) == pytest.approx(
        probabilistic_sharpe_ratio(rets, 0.0), abs=1e-12)


def test_dsr_deflates_for_multiple_trials():
    rets = _sample()
    assert deflated_sharpe_ratio(rets, 10) < probabilistic_sharpe_ratio(rets, 0.0)


def test_dsr_decreasing_in_trials():
    rets = _sample()
    assert deflated_sharpe_ratio(rets, 100) < deflated_sharpe_ratio(rets, 10) \
        < deflated_sharpe_ratio(rets, 2)


def test_dsr_in_unit_interval():
    rets = _sample()
    assert all(0 <= deflated_sharpe_ratio(rets, n) <= 1 for n in (1, 5, 50, 500))


def test_dsr_validation():
    with pytest.raises(ValueError):
        deflated_sharpe_ratio(_sample(), 0)


def test_validation():
    with pytest.raises(ValueError):
        probabilistic_sharpe_ratio([0.01] * 10)   # zero variance
    with pytest.raises(ValueError):
        probabilistic_sharpe_ratio([0.01])         # too few
