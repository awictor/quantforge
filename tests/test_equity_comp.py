"""Warrants and employee stock options."""

import math

import pytest

from quantforge import (
    dilution_factor, warrant_price, eso_expected_life, eso_value,
    pv_dividends, discrete_dividend_price, forward_with_dividends,
    conversion_value, straight_bond_floor, convertible_bond_value,
)
from quantforge.bsm import call_price


DIVS = [(0.25, 2.0), (0.75, 2.0)]


S, K, T, R, SIG = 50.0, 50.0, 5.0, 0.05, 0.3


def test_dilution_factor():
    assert dilution_factor(1_000_000, 100_000) == pytest.approx(1e6 / 1.1e6)
    assert dilution_factor(1e6, 0) == 1.0


def test_warrant_below_vanilla_call():
    van = call_price(S, K, T, R, SIG)
    w = warrant_price(S, K, T, R, SIG, 1_000_000, 100_000)
    assert w < van
    assert w / van == pytest.approx(dilution_factor(1_000_000, 100_000))


def test_expected_life_between_vesting_and_term():
    life = eso_expected_life(2.0, 10.0, 0.15)
    assert 2.0 < life < 10.0


def test_expected_life_no_exit_is_full_term():
    assert eso_expected_life(2.0, 10.0, 0.0) == 10.0


def test_eso_below_vanilla_on_full_term():
    eso = eso_value(S, K, 10.0, R, SIG, 2.0, 0.15, 0.03)
    assert eso < call_price(S, K, 10.0, R, SIG)


def test_eso_reduces_to_vanilla():
    assert eso_value(S, K, 10.0, R, SIG, 2.0, 0.0, 0.0) == pytest.approx(
        call_price(S, K, 10.0, R, SIG), abs=1e-9)


def test_higher_exit_rate_lowers_value():
    assert eso_value(S, K, 10.0, R, SIG, 2.0, 0.3, 0.0) < \
        eso_value(S, K, 10.0, R, SIG, 2.0, 0.1, 0.0)


def test_forfeiture_lowers_value():
    assert eso_value(S, K, 10.0, R, SIG, 2.0, 0.15, 0.1) < \
        eso_value(S, K, 10.0, R, SIG, 2.0, 0.15, 0.0)


def test_conversion_value():
    assert conversion_value(50, 20) == 1000


def test_bond_floor_below_face():
    floor = straight_bond_floor(1000, 0.04, 5, 0.05, 0.01)
    assert floor < 1000


def test_convertible_above_floor_and_parity():
    floor = straight_bond_floor(1000, 0.04, 5, 0.05, 0.3 * 0 + 0.01)
    c = convertible_bond_value(50, 20, 1000, 0.04, 5, 0.05, 0.3, 0.01)
    assert c >= floor
    assert c >= conversion_value(50, 20) - 1e-6


def test_convertible_deep_itm_approaches_parity():
    c = convertible_bond_value(200, 20, 1000, 0.04, 5, 0.05, 0.3, 0.01)
    conv = conversion_value(200, 20)
    assert abs(c - conv) / conv < 0.15


def test_convertible_deep_otm_approaches_floor():
    floor = straight_bond_floor(1000, 0.04, 5, 0.05, 0.01)
    c = convertible_bond_value(5, 20, 1000, 0.04, 5, 0.05, 0.3, 0.01)
    assert abs(c - floor) / floor < 0.05


def test_convertible_spread_and_vol_effects():
    base = convertible_bond_value(50, 20, 1000, 0.04, 5, 0.05, 0.3, 0.01)
    assert convertible_bond_value(50, 20, 1000, 0.04, 5, 0.05, 0.3, 0.05) < base
    assert convertible_bond_value(50, 20, 1000, 0.04, 5, 0.05, 0.5, 0.01) > base


def test_convertible_validation():
    with pytest.raises(ValueError):
        convertible_bond_value(50, 0, 1000, 0.04, 5, 0.05, 0.3)


def test_pv_dividends():
    assert pv_dividends(DIVS, 0.05) == pytest.approx(
        2 * math.exp(-0.05 * 0.25) + 2 * math.exp(-0.05 * 0.75))


def test_no_dividends_is_vanilla():
    assert discrete_dividend_price(100, 100, 1.0, 0.05, 0.25, []) == pytest.approx(
        call_price(100, 100, 1.0, 0.05, 0.25), abs=1e-9)


def test_dividends_lower_call_and_match_adjusted_spot():
    c = discrete_dividend_price(100, 100, 1.0, 0.05, 0.25, DIVS)
    assert c < call_price(100, 100, 1.0, 0.05, 0.25)
    pv = pv_dividends(DIVS, 0.05)
    assert c == pytest.approx(call_price(100 - pv, 100, 1.0, 0.05, 0.25, b=0.05),
                              abs=1e-12)


def test_discrete_dividend_parity():
    c = discrete_dividend_price(100, 100, 1.0, 0.05, 0.25, DIVS, True)
    p = discrete_dividend_price(100, 100, 1.0, 0.05, 0.25, DIVS, False)
    pv = pv_dividends(DIVS, 0.05)
    assert c - p == pytest.approx((100 - pv) - 100 * math.exp(-0.05), abs=1e-9)


def test_forward_with_dividends():
    pv = pv_dividends(DIVS, 0.05)
    assert forward_with_dividends(100, 1.0, 0.05, DIVS) == pytest.approx(
        (100 - pv) * math.exp(0.05))
    assert forward_with_dividends(100, 1.0, 0.05, DIVS) < 100 * math.exp(0.05)


def test_dividends_after_expiry_excluded():
    assert discrete_dividend_price(100, 100, 1.0, 0.05, 0.25, [(2.0, 5.0)]) == \
        pytest.approx(call_price(100, 100, 1.0, 0.05, 0.25), abs=1e-9)


def test_dividend_validation():
    with pytest.raises(ValueError):
        discrete_dividend_price(1, 100, 1.0, 0.05, 0.25, [(0.5, 100)])


def test_validation():
    with pytest.raises(ValueError):
        dilution_factor(-1, 100)
    with pytest.raises(ValueError):
        eso_expected_life(11, 10, 0.1)
    with pytest.raises(ValueError):
        eso_value(S, K, 10.0, R, SIG, 2.0, 0.15, -1)
