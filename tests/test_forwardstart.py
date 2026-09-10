"""Tests for forward-start and cliquet options."""

import math

import pytest

from quantforge import (
    forward_start_price, cliquet_price, call_price, put_price, OptionType,
)


def test_forward_start_reduces_to_vanilla_at_tstart_zero():
    # As t_start -> 0 the strike is set at alpha*S today, i.e. a vanilla struck
    # at alpha*S with the full time to expiry.
    S, T, r, sigma, alpha = 100.0, 1.0, 0.05, 0.25, 1.0
    fs = forward_start_price(S, 1e-9, T, r, sigma, alpha, OptionType.CALL)
    vanilla = call_price(S, alpha * S, T, r, sigma)
    assert fs == pytest.approx(vanilla, abs=1e-3)


def test_forward_start_scales_linearly_in_spot():
    # BSM homogeneity: doubling spot doubles a forward-start value.
    kw = dict(t_start=0.25, t_expiry=1.0, r=0.05, sigma=0.3, alpha=1.0)
    v1 = forward_start_price(100.0, **kw)
    v2 = forward_start_price(200.0, **kw)
    assert v2 == pytest.approx(2.0 * v1, rel=1e-9)


def test_forward_start_independent_of_strike_level():
    # The value depends on alpha (moneyness) but not on the absolute future
    # strike, so two spots give values that scale exactly by the spot ratio.
    kw = dict(t_start=0.5, t_expiry=1.5, r=0.03, sigma=0.2, alpha=1.1,
              option_type=OptionType.PUT)
    v_a = forward_start_price(80.0, **kw)
    v_b = forward_start_price(130.0, **kw)
    assert v_b / v_a == pytest.approx(130.0 / 80.0, rel=1e-9)


def test_atm_forward_start_positive():
    v = forward_start_price(100, 0.5, 1.0, 0.05, 0.25, alpha=1.0)
    assert v > 0


def test_cliquet_equals_sum_of_periods():
    S, r, sigma = 100.0, 0.05, 0.2
    resets = [0.25, 0.5, 0.75, 1.0]
    total = cliquet_price(S, resets, r, sigma, alpha=1.0, option_type=OptionType.CALL)

    # Reconstruct: first period is a vanilla struck at S; the rest are
    # forward-starts.
    manual = call_price(S, S, 0.25, r, sigma)
    starts = [0.25, 0.5, 0.75]
    ends = [0.5, 0.75, 1.0]
    for ts, te in zip(starts, ends):
        manual += forward_start_price(S, ts, te, r, sigma, 1.0, OptionType.CALL)
    assert total == pytest.approx(manual, rel=1e-12)


def test_cliquet_more_periods_costs_more():
    S, r, sigma = 100.0, 0.05, 0.2
    two = cliquet_price(S, [0.5, 1.0], r, sigma)
    four = cliquet_price(S, [0.25, 0.5, 0.75, 1.0], r, sigma)
    assert four > two  # more resets = more optionality


def test_invalid_inputs_rejected():
    with pytest.raises(ValueError):
        forward_start_price(100, 1.0, 0.5, 0.05, 0.2)  # t_start >= t_expiry
    with pytest.raises(ValueError):
        forward_start_price(100, 0.25, 1.0, 0.05, 0.2, alpha=-1)  # bad alpha
    with pytest.raises(ValueError):
        cliquet_price(100, [0.5, 0.5], 0.05, 0.2)  # not strictly increasing
    with pytest.raises(ValueError):
        cliquet_price(100, [], 0.05, 0.2)  # empty schedule
