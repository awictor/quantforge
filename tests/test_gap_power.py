"""Tests for gap and power options."""

import math

import pytest

from quantforge import gap_option, power_option, call_price, put_price, OptionType


# --- Gap options ---
def test_gap_equal_strikes_is_vanilla_call():
    v = gap_option(100, 100, 100, 1.0, 0.05, 0.2, OptionType.CALL)
    assert v == pytest.approx(call_price(100, 100, 1.0, 0.05, 0.2), abs=1e-9)


def test_gap_equal_strikes_is_vanilla_put():
    v = gap_option(100, 100, 100, 1.0, 0.05, 0.2, OptionType.PUT)
    assert v == pytest.approx(put_price(100, 100, 1.0, 0.05, 0.2), abs=1e-9)


def test_gap_payoff_can_be_negative_valued():
    # Trigger below payoff strike: pays S - K_payoff whenever S > K_trigger, so
    # the region K_trigger < S < K_payoff contributes negative payoff -> the
    # option is worth less than the vanilla struck at K_trigger.
    gap = gap_option(100, 90, 110, 1.0, 0.05, 0.25, OptionType.CALL)
    vanilla = call_price(100, 90, 1.0, 0.05, 0.25)
    assert gap < vanilla


def test_gap_higher_payoff_strike_lowers_call():
    lo = gap_option(100, 100, 100, 1.0, 0.05, 0.25, OptionType.CALL)
    hi = gap_option(100, 100, 120, 1.0, 0.05, 0.25, OptionType.CALL)
    assert hi < lo


# --- Power options ---
def test_power_one_is_vanilla():
    v = power_option(100, 100, 1.0, 0.05, 0.2, power=1.0, option_type=OptionType.CALL)
    assert v == pytest.approx(call_price(100, 100, 1.0, 0.05, 0.2), abs=1e-9)


def test_power_two_matches_monte_carlo():
    # Cross-checked against a terminal-value Monte Carlo (~2433).
    v = power_option(100, 10000, 1.0, 0.05, 0.2, power=2.0, option_type=OptionType.CALL)
    assert v == pytest.approx(2433, rel=0.02)


def test_power_put_call_relationship_positive():
    c = power_option(100, 10000, 1.0, 0.05, 0.2, power=2.0, option_type=OptionType.CALL)
    p = power_option(100, 10000, 1.0, 0.05, 0.2, power=2.0, option_type=OptionType.PUT)
    assert c > 0 and p > 0


def test_power_zero_time_intrinsic():
    v = power_option(100, 5000, 0.0, 0.05, 0.2, power=2.0, option_type=OptionType.CALL)
    assert v == pytest.approx(max(100 ** 2 - 5000, 0.0), abs=1e-6)


def test_rejects_bad_inputs():
    with pytest.raises(ValueError):
        power_option(100, 100, 1.0, 0.05, 0.2, power=-1.0)
    with pytest.raises(ValueError):
        gap_option(100, 100, -5, 1.0, 0.05, 0.2)
