"""Tests for digital super-replication by a vanilla spread."""

import pytest

from quantforge import (
    digital_call_overhedge, digital_put_overhedge, overhedge_payoff,
    cash_or_nothing, OptionType,
)


def test_call_overhedge_cost_exceeds_digital():
    # The dominating spread must cost at least the fair digital value.
    oh = digital_call_overhedge(100, 100, 0.25, 0.05, 0.2, cash=1.0, width=2)
    assert oh.cushion > 0
    assert oh.cost > oh.digital_value


def test_put_overhedge_cost_exceeds_digital():
    oh = digital_put_overhedge(100, 100, 0.25, 0.05, 0.2, cash=1.0, width=2)
    assert oh.cushion > 0


def test_cushion_shrinks_with_width():
    wide = digital_call_overhedge(100, 100, 0.25, 0.05, 0.2, width=10)
    narrow = digital_call_overhedge(100, 100, 0.25, 0.05, 0.2, width=0.5)
    assert wide.cushion > narrow.cushion


def test_cost_converges_to_digital_as_width_shrinks():
    oh = digital_call_overhedge(100, 100, 0.25, 0.05, 0.2, width=0.05)
    assert oh.cost == pytest.approx(oh.digital_value, abs=1e-3)


@pytest.mark.parametrize("ST", [80, 95, 99, 99.999, 100.001, 101, 110, 130])
def test_call_spread_dominates_digital_payoff(ST):
    cash = 1.0
    oh = digital_call_overhedge(100, 100, 0.25, 0.05, 0.2, cash=cash, width=2)
    pay = overhedge_payoff(oh, ST, is_call=True)
    digital = cash if ST > 100 else 0.0
    assert pay >= digital - 1e-9


@pytest.mark.parametrize("ST", [70, 90, 99, 100.001, 105, 120])
def test_put_spread_dominates_digital_payoff(ST):
    cash = 1.0
    oh = digital_put_overhedge(100, 100, 0.25, 0.05, 0.2, cash=cash, width=2)
    pay = overhedge_payoff(oh, ST, is_call=False)
    digital = cash if ST < 100 else 0.0
    assert pay >= digital - 1e-9


def test_call_spread_payoff_caps_at_cash():
    oh = digital_call_overhedge(100, 100, 0.25, 0.05, 0.2, cash=1.0, width=2)
    assert overhedge_payoff(oh, 200, is_call=True) == pytest.approx(1.0)
    assert overhedge_payoff(oh, 50, is_call=True) == pytest.approx(0.0)


def test_digital_value_matches_cash_or_nothing():
    oh = digital_call_overhedge(100, 105, 0.5, 0.04, 0.3, cash=2.0, width=1)
    ref = cash_or_nothing(100, 105, 0.5, 0.04, 0.3, OptionType.CALL, cash=2.0)
    assert oh.digital_value == pytest.approx(ref, abs=1e-12)


def test_rejects_bad_width():
    with pytest.raises(ValueError):
        digital_call_overhedge(100, 100, 0.25, 0.05, 0.2, width=-1)
