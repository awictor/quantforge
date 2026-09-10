"""Tests for smile-aware delta (sticky-strike vs sticky-delta)."""

import math

import pytest

from quantforge import (
    delta, vega, smile_delta, smile_delta_from_smile, skew_slope,
    StickyRule, OptionType,
)


BASE = dict(S=100.0, K=105.0, t=0.5, r=0.04, sigma=0.25)


def test_sticky_strike_equals_bs_delta():
    d = smile_delta(**BASE, option_type=OptionType.CALL, sticky=StickyRule.STRIKE)
    assert d == pytest.approx(delta(**BASE, option_type=OptionType.CALL), abs=1e-12)


def test_sticky_delta_flat_smile_equals_bs():
    # Zero skew slope -> no adjustment -> equals BS delta.
    d = smile_delta(**BASE, dsigma_dk=0.0, option_type=OptionType.CALL,
                    sticky=StickyRule.DELTA)
    assert d == pytest.approx(delta(**BASE, option_type=OptionType.CALL), abs=1e-12)


def test_sticky_delta_adds_vega_skew_term():
    slope = -0.5  # downward (equity) skew: vol falls as k rises
    d_bs = delta(**BASE, option_type=OptionType.CALL)
    v = vega(**BASE)
    d = smile_delta(**BASE, dsigma_dk=slope, option_type=OptionType.CALL,
                    sticky=StickyRule.DELTA)
    assert d == pytest.approx(d_bs - v * slope / BASE["S"], abs=1e-12)


def test_negative_skew_raises_call_delta():
    # With a downward skew (dsigma/dk < 0) the sticky-delta call delta exceeds
    # the BS delta: -vega*(negative)/S > 0.
    d_bs = delta(**BASE, option_type=OptionType.CALL)
    d = smile_delta(**BASE, dsigma_dk=-0.4, option_type=OptionType.CALL,
                    sticky=StickyRule.DELTA)
    assert d > d_bs


def test_skew_slope_of_linear_smile():
    # A smile linear in k = ln(K/F): sigma(K) = 0.2 - 0.3 * ln(K/F).
    F = 100.0
    smile = lambda K: 0.2 - 0.3 * math.log(K / F)
    slope = skew_slope(smile, K=110.0, F=F)
    assert slope == pytest.approx(-0.3, abs=1e-6)


def test_smile_delta_from_smile_matches_manual():
    F = 100.0 * math.exp(0.04 * 0.5)
    smile = lambda K: 0.25 - 0.2 * math.log(K / F)
    K = 105.0
    d_auto = smile_delta_from_smile(100.0, K, 0.5, 0.04, smile,
                                    option_type=OptionType.CALL)
    # Manual: same sigma and slope fed into smile_delta.
    sigma = smile(K)
    slope = skew_slope(smile, K, F)
    d_manual = smile_delta(100.0, K, 0.5, 0.04, sigma, dsigma_dk=slope,
                           option_type=OptionType.CALL, sticky=StickyRule.DELTA)
    assert d_auto == pytest.approx(d_manual, abs=1e-12)


def test_put_sticky_delta_adjustment_sign():
    slope = -0.4
    d_bs = delta(**BASE, option_type=OptionType.PUT)
    v = vega(**BASE)
    d = smile_delta(**BASE, dsigma_dk=slope, option_type=OptionType.PUT,
                    sticky=StickyRule.DELTA)
    assert d == pytest.approx(d_bs - v * slope / BASE["S"], abs=1e-12)
