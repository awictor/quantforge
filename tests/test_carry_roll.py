"""Tests for carry-roll / roll-down P&L."""

import math

import pytest

from quantforge import carry_roll_pnl, price, OptionType


def test_long_option_bleeds_over_the_roll():
    r = carry_roll_pnl(100, 100, 1.0, 0.05, 0.2, horizon=0.25, option_type=OptionType.CALL)
    assert r.theta_roll < 0                 # a long ATM call loses value rolling forward
    assert r.value_rolled_static < r.value_now


def test_short_position_flips_sign():
    lo = carry_roll_pnl(100, 100, 1.0, 0.05, 0.2, horizon=0.25, qty=1)
    sh = carry_roll_pnl(100, 100, 1.0, 0.05, 0.2, horizon=0.25, qty=-1)
    assert sh.theta_roll == pytest.approx(-lo.theta_roll)
    assert sh.theta_roll > 0


def test_forward_spot_is_carry_grown():
    r = carry_roll_pnl(100, 100, 1.0, 0.05, 0.2, horizon=0.5, b=0.03)
    assert r.forward_spot == pytest.approx(100 * math.exp(0.03 * 0.5))


def test_rolled_value_matches_reprice():
    r = carry_roll_pnl(100, 105, 1.0, 0.04, 0.25, horizon=0.4, option_type=OptionType.PUT)
    expected = price(r.forward_spot, 105, 1.0 - 0.4, 0.04, 0.25, OptionType.PUT, b=0.04)
    assert r.value_rolled_static == pytest.approx(expected, abs=1e-9)


def test_value_now_matches_price():
    r = carry_roll_pnl(100, 100, 1.0, 0.05, 0.2, horizon=0.25)
    assert r.value_now == pytest.approx(price(100, 100, 1.0, 0.05, 0.2, OptionType.CALL))


def test_scales_with_quantity():
    one = carry_roll_pnl(100, 100, 1.0, 0.05, 0.2, horizon=0.25, qty=1)
    ten = carry_roll_pnl(100, 100, 1.0, 0.05, 0.2, horizon=0.25, qty=10)
    assert ten.theta_roll == pytest.approx(10 * one.theta_roll)


def test_rejects_bad_horizon():
    with pytest.raises(ValueError):
        carry_roll_pnl(100, 100, 1.0, 0.05, 0.2, horizon=0.0)
    with pytest.raises(ValueError):
        carry_roll_pnl(100, 100, 1.0, 0.05, 0.2, horizon=1.5)
