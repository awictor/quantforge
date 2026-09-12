"""Probabilistic Sharpe ratio and minimum track record length."""

import random
import statistics

import pytest

from quantforge import probabilistic_sharpe_ratio, minimum_track_record_length


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


def test_validation():
    with pytest.raises(ValueError):
        probabilistic_sharpe_ratio([0.01] * 10)   # zero variance
    with pytest.raises(ValueError):
        probabilistic_sharpe_ratio([0.01])         # too few
