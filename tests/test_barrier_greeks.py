"""Tests for single-barrier option Greeks (finite differences)."""

import pytest

from quantforge import (
    barrier_greeks, barrier_option, greeks, Barrier, OptionType,
)


def test_far_down_out_call_greeks_match_vanilla():
    # A barrier far below spot is (almost) never hit, so the Greeks match the
    # vanilla call's Greeks.
    g = barrier_greeks(100, 100, 1.0, 1.0, 0.05, 0.2, OptionType.CALL,
                       Barrier.DOWN_OUT)
    e = greeks(100, 100, 1.0, 0.05, 0.2, OptionType.CALL)
    assert g["delta"] == pytest.approx(e.delta, abs=1e-4)
    assert g["gamma"] == pytest.approx(e.gamma, abs=1e-4)
    assert g["vega"] == pytest.approx(e.vega, abs=1e-2)
    assert g["theta"] == pytest.approx(e.theta, abs=1e-2)


def test_price_field_matches_direct():
    g = barrier_greeks(100, 90, 95, 0.5, 0.08, 0.25, OptionType.CALL,
                       Barrier.DOWN_OUT, b=0.04, rebate=3.0)
    assert g["price"] == pytest.approx(
        barrier_option(100, 90, 95, 0.5, 0.08, 0.25, OptionType.CALL,
                       Barrier.DOWN_OUT, b=0.04, rebate=3.0), abs=1e-9)


def test_delta_matches_re_difference():
    S, K, H = 100, 100, 90
    g = barrier_greeks(S, K, H, 0.5, 0.05, 0.25, OptionType.CALL, Barrier.DOWN_OUT)
    h = 0.05
    up = barrier_option(S + h, K, H, 0.5, 0.05, 0.25, OptionType.CALL, Barrier.DOWN_OUT)
    dn = barrier_option(S - h, K, H, 0.5, 0.05, 0.25, OptionType.CALL, Barrier.DOWN_OUT)
    assert g["delta"] == pytest.approx((up - dn) / (2 * h), abs=1e-2)


def test_all_fields_present():
    g = barrier_greeks(100, 100, 110, 0.5, 0.05, 0.2, OptionType.CALL, Barrier.UP_OUT)
    assert set(g) == {"price", "delta", "gamma", "vega", "theta"}


def test_knock_in_plus_knock_out_delta_equals_vanilla():
    # Delta is linear in the price, so in-out parity carries to delta.
    S, K, H = 100, 100, 90
    ki = barrier_greeks(S, K, H, 1.0, 0.05, 0.25, OptionType.CALL, Barrier.DOWN_IN)
    ko = barrier_greeks(S, K, H, 1.0, 0.05, 0.25, OptionType.CALL, Barrier.DOWN_OUT)
    van = greeks(S, K, 1.0, 0.05, 0.25, OptionType.CALL)
    assert ki["delta"] + ko["delta"] == pytest.approx(van.delta, abs=1e-3)
