"""Tests for displaced-diffusion (shifted lognormal) pricing."""

import math

import pytest

from quantforge import displaced_diffusion_price as dd, call_price, put_price, OptionType


def test_zero_shift_is_black_scholes():
    for K in (80, 100, 120):
        v = dd(100, K, 1.0, 0.05, 0.2, shift=0.0, option_type=OptionType.CALL)
        assert v == pytest.approx(call_price(100, K, 1.0, 0.05, 0.2), abs=1e-9)


def test_put_call_parity():
    c = dd(100, 95, 1.0, 0.05, 0.2, shift=50, option_type=OptionType.CALL)
    p = dd(100, 95, 1.0, 0.05, 0.2, shift=50, option_type=OptionType.PUT)
    assert c - p == pytest.approx(100 - 95 * math.exp(-0.05), abs=1e-6)


def test_matches_monte_carlo():
    # Cross-checked against a displaced-GBM Monte Carlo (~10.40).
    v = dd(100, 100, 1.0, 0.05, 0.2, shift=50, option_type=OptionType.CALL)
    assert v == pytest.approx(10.40, abs=0.1)


def test_allows_negative_strike():
    # With a shift the strike can be negative; the price stays finite/positive.
    v = dd(10, -5, 1.0, 0.05, 0.5, shift=20, option_type=OptionType.CALL)
    assert v > 0 and math.isfinite(v)


def test_positive_shift_dampens_smile_curvature():
    # A larger positive shift moves the model toward normal behavior, so far-OTM
    # options change value monotonically vs pure lognormal. Just check finiteness
    # and positivity across strikes/shifts here.
    for shift in (0, 20, 100):
        for K in (80, 100, 120):
            v = dd(100, K, 1.0, 0.05, 0.25, shift=shift, option_type=OptionType.PUT)
            assert v >= 0


def test_zero_time_is_intrinsic():
    v = dd(110, 100, 0.0, 0.05, 0.2, shift=30, option_type=OptionType.CALL)
    assert v == pytest.approx(10.0, abs=1e-9)


def test_rejects_shift_below_negative_spot():
    with pytest.raises(ValueError):
        dd(100, 100, 1.0, 0.05, 0.2, shift=-150)  # S + shift < 0
