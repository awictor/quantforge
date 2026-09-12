"""Callable and puttable bond pricing on a short-rate tree."""

import pytest

from quantforge import (
    callable_bond_price, straight_bond_tree_price, call_option_value,
    callable_bond_price_with_spread, option_adjusted_spread,
)


FACE, CPN, MAT, R0, SIG = 100.0, 0.05, 5, 0.04, 0.15


def test_callable_below_straight():
    straight = straight_bond_tree_price(FACE, CPN, MAT, R0, SIG)
    callable_ = callable_bond_price(FACE, CPN, MAT, R0, SIG, call_price=100)
    assert callable_ <= straight + 1e-9


def test_puttable_above_straight():
    straight = straight_bond_tree_price(FACE, CPN, MAT, R0, SIG)
    puttable = callable_bond_price(FACE, CPN, MAT, R0, SIG, put_price=100)
    assert puttable >= straight - 1e-9


def test_no_optionality_is_straight():
    assert callable_bond_price(FACE, CPN, MAT, R0, SIG) == pytest.approx(
        straight_bond_tree_price(FACE, CPN, MAT, R0, SIG), abs=1e-12)


def test_call_option_value_non_negative():
    assert call_option_value(FACE, CPN, MAT, R0, SIG, 100) >= 0


def test_higher_vol_lowers_callable():
    lo = callable_bond_price(FACE, CPN, MAT, R0, 0.15, call_price=100)
    hi = callable_bond_price(FACE, CPN, MAT, R0, 0.30, call_price=100)
    assert hi < lo


def test_higher_call_price_raises_callable():
    base = callable_bond_price(FACE, CPN, MAT, R0, SIG, call_price=100)
    higher = callable_bond_price(FACE, CPN, MAT, R0, SIG, call_price=110)
    assert higher > base


def test_spread_zero_matches_base():
    base = callable_bond_price(FACE, CPN, MAT, R0, SIG, call_price=110)
    assert callable_bond_price_with_spread(FACE, CPN, MAT, R0, SIG, 0.0,
                                           call_price=110) == pytest.approx(base, abs=1e-12)


def test_positive_spread_lowers_price():
    base = callable_bond_price_with_spread(FACE, CPN, MAT, R0, SIG, 0.0, call_price=110)
    assert callable_bond_price_with_spread(FACE, CPN, MAT, R0, SIG, 0.01,
                                           call_price=110) < base


def test_oas_recovers_spread():
    target = callable_bond_price_with_spread(FACE, CPN, MAT, R0, SIG, 0.02, call_price=110)
    assert option_adjusted_spread(target, FACE, CPN, MAT, R0, SIG,
                                  call_price=110) == pytest.approx(0.02, abs=1e-6)


def test_oas_zero_at_model_price():
    base = callable_bond_price(FACE, CPN, MAT, R0, SIG, call_price=110)
    assert option_adjusted_spread(base, FACE, CPN, MAT, R0, SIG,
                                  call_price=110) == pytest.approx(0.0, abs=1e-6)


def test_oas_positive_below_model():
    base = callable_bond_price(FACE, CPN, MAT, R0, SIG, call_price=110)
    assert option_adjusted_spread(base * 0.95, FACE, CPN, MAT, R0, SIG,
                                  call_price=110) > 0


def test_oas_validation():
    with pytest.raises(ValueError):
        option_adjusted_spread(-5, FACE, CPN, MAT, R0, SIG)


def test_validation():
    with pytest.raises(ValueError):
        callable_bond_price(-1, CPN, MAT, R0, SIG)
    with pytest.raises(ValueError):
        callable_bond_price(FACE, CPN, 0, R0, SIG)
