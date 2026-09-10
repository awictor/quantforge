"""Tests for standalone barrier rebate cashflows."""

import math

import pytest

from quantforge import barrier_rebate, one_touch, no_touch


def test_knockout_rebate_is_one_touch():
    v = barrier_rebate(100, 120, 1.0, 0.05, 0.2, knock="out", payoff_at_hit=True)
    assert v == pytest.approx(one_touch(100, 120, 1.0, 0.05, 0.2, cash=1.0))


def test_knockin_rebate_is_no_touch():
    v = barrier_rebate(100, 120, 1.0, 0.05, 0.2, knock="in")
    assert v == pytest.approx(no_touch(100, 120, 1.0, 0.05, 0.2, cash=1.0))


def test_pay_at_hit_at_least_pay_at_expiry():
    hit = barrier_rebate(100, 120, 1.0, 0.05, 0.2, knock="out", payoff_at_hit=True)
    exp = barrier_rebate(100, 120, 1.0, 0.05, 0.2, knock="out", payoff_at_hit=False)
    assert hit >= exp


def test_knockout_expiry_plus_knockin_is_discounted_cash():
    ko = barrier_rebate(100, 120, 1.0, 0.05, 0.2, knock="out", payoff_at_hit=False)
    ki = barrier_rebate(100, 120, 1.0, 0.05, 0.2, knock="in")
    assert ko + ki == pytest.approx(math.exp(-0.05), abs=1e-9)


def test_cash_scales():
    a = barrier_rebate(100, 120, 1.0, 0.05, 0.2, knock="out", cash=1.0)
    b = barrier_rebate(100, 120, 1.0, 0.05, 0.2, knock="out", cash=25.0)
    assert b == pytest.approx(25.0 * a)


def test_rejects_bad_knock():
    with pytest.raises(ValueError):
        barrier_rebate(100, 120, 1.0, 0.05, 0.2, knock="maybe")
