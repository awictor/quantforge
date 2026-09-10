"""Tests for one-touch / no-touch binary options."""

import math

import pytest

from quantforge import one_touch, no_touch, cash_or_nothing, OptionType


def test_one_touch_plus_no_touch_is_discounted_cash():
    # With payment at expiry, touch and no-touch partition the cash.
    S, t, r, sigma = 100, 1.0, 0.05, 0.2
    for H in (120, 80):
        ot = one_touch(S, H, t, r, sigma, payoff_at_hit=False)
        nt = no_touch(S, H, t, r, sigma)
        assert ot + nt == pytest.approx(math.exp(-r * t), abs=1e-9)


def test_pay_at_hit_at_least_pay_at_expiry():
    # Receiving the cash earlier (on hit) is worth more than at expiry.
    S, t, r, sigma = 100, 1.0, 0.05, 0.2
    hit_now = one_touch(S, 120, t, r, sigma, payoff_at_hit=True)
    hit_exp = one_touch(S, 120, t, r, sigma, payoff_at_hit=False)
    assert hit_now >= hit_exp


def test_touch_probability_in_unit_interval():
    S, t, r, sigma = 100, 1.0, 0.05, 0.3
    for H in (70, 90, 110, 140):
        p = one_touch(S, H, t, r, sigma, payoff_at_hit=False) * math.exp(r * t)
        assert 0.0 <= p <= 1.0


def test_closer_barrier_more_likely_touched():
    S, t, r, sigma = 100, 1.0, 0.05, 0.2
    near = one_touch(S, 105, t, r, sigma, payoff_at_hit=False)
    far = one_touch(S, 140, t, r, sigma, payoff_at_hit=False)
    assert near > far


def test_reference_value_matches_monte_carlo():
    # Cross-checked against a fine barrier-crossing Monte Carlo (~0.40).
    v = one_touch(100, 120, 1.0, 0.05, 0.2, payoff_at_hit=True)
    assert v == pytest.approx(0.40, abs=0.02)


def test_cash_scales_linearly():
    S, t, r, sigma, H = 100, 1.0, 0.05, 0.2, 120
    unit = one_touch(S, H, t, r, sigma, cash=1.0)
    ten = one_touch(S, H, t, r, sigma, cash=10.0)
    assert ten == pytest.approx(10.0 * unit)


def test_no_touch_relates_to_double_of_cash_or_nothing():
    # A no-touch that never gets near the barrier approaches discounted cash.
    S, t, r, sigma = 100, 1.0, 0.05, 0.2
    nt = no_touch(S, 1000.0, t, r, sigma)  # far away barrier
    assert nt == pytest.approx(math.exp(-r * t), abs=1e-3)


def test_zero_time_no_touch_pays_full_cash():
    assert no_touch(100, 120, 0.0, 0.05, 0.2) == pytest.approx(1.0)
    assert one_touch(100, 120, 0.0, 0.05, 0.2) == pytest.approx(0.0)
