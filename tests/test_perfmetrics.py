"""Performance and drawdown statistics (perfmetrics module)."""

import math
import statistics

import pytest

from quantforge import (
    sharpe_ratio, sortino_ratio, max_drawdown, calmar_ratio, hit_rate,
    profit_factor,
)


R = [0.01, -0.005, 0.02, 0.008, -0.012, 0.015, 0.003, -0.002]


def test_sharpe_matches_manual():
    m = sum(R) / len(R)
    sd = statistics.stdev(R)
    assert sharpe_ratio(R) == pytest.approx(m / sd * math.sqrt(252), abs=1e-9)


def test_sortino_above_sharpe():
    assert sortino_ratio(R) > sharpe_ratio(R)


def test_max_drawdown_known_series():
    # peak 1.1 then trough 0.55 -> 50% drawdown.
    assert max_drawdown([0.1, -0.5, 0.2]) == pytest.approx(0.5, abs=1e-12)


def test_max_drawdown_all_positive_is_zero():
    assert max_drawdown([0.01, 0.02, 0.03]) == 0.0


def test_hit_rate():
    assert hit_rate(R) == pytest.approx(5 / 8, abs=1e-12)


def test_profit_factor_and_infinite_case():
    assert profit_factor(R) > 1.0
    assert profit_factor([0.01, 0.02]) == float("inf")


def test_calmar_positive_and_no_drawdown_raises():
    assert calmar_ratio(R) > 0.0
    with pytest.raises(ValueError):
        calmar_ratio([0.01, 0.02])   # no drawdown


def test_validation():
    with pytest.raises(ValueError):
        sharpe_ratio([0.01])          # need >= 2
    with pytest.raises(ValueError):
        sharpe_ratio([0.01, 0.01])    # zero variance
    with pytest.raises(ValueError):
        sortino_ratio([0.01, 0.02])   # no downside
    with pytest.raises(ValueError):
        hit_rate([])
