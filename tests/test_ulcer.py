"""Ulcer index and pain-based drawdown ratios."""

import pytest

from quantforge import (
    ulcer_index, pain_index, ulcer_performance_index, pain_ratio, max_drawdown,
)


RETS = [0.05, -0.10, 0.02, -0.08, 0.03, 0.06, -0.04, 0.01, -0.12, 0.05]


def test_no_drawdown_zero():
    up = [0.01] * 10
    assert ulcer_index(up) == 0.0
    assert pain_index(up) == 0.0


def test_ulcer_at_least_pain():
    # RMS >= mean for the non-negative drawdown series.
    assert ulcer_index(RETS) >= pain_index(RETS)


def test_bounded_by_max_drawdown():
    mdd = max_drawdown(RETS)
    assert ulcer_index(RETS) <= mdd + 1e-9
    assert pain_index(RETS) <= mdd + 1e-9


def test_ratios_finite_with_drawdown():
    assert ulcer_performance_index(RETS) == pytest.approx(
        (sum(RETS) / len(RETS) * 252) / ulcer_index(RETS))
    assert pain_ratio(RETS) == pytest.approx(
        (sum(RETS) / len(RETS) * 252) / pain_index(RETS))


def test_no_drawdown_ratios_raise():
    up = [0.01] * 10
    with pytest.raises(ValueError):
        ulcer_performance_index(up)
    with pytest.raises(ValueError):
        pain_ratio(up)
