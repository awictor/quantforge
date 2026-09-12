"""Par yields and par swap rates from a discount curve."""

import math

import pytest

from quantforge import par_yield, par_bond_price


def _flat(r):
    return lambda t: (1 + r) ** (-t)


def test_flat_curve_par_equals_flat_rate():
    P = _flat(0.04)
    for T in (1, 5, 10, 30):
        assert abs(par_yield(P, T) - 0.04) < 1e-9


def test_par_bond_prices_to_face():
    P = _flat(0.04)
    y = par_yield(P, 5)
    assert abs(par_bond_price(P, y, 5) - 100.0) < 1e-9


def test_one_year_par_equals_spot():
    P = _flat(0.04)
    assert abs(par_yield(P, 1) - 0.04) < 1e-9


def test_semiannual_par_bond_to_face():
    P = _flat(0.04)
    y = par_yield(P, 5, freq=2)
    assert abs(par_bond_price(P, y, 5, freq=2) - 100.0) < 1e-9


def test_upward_curve_par_between_zeros():
    Pu = lambda t: math.exp(-(0.02 + 0.003 * t) * t) if t > 0 else 1.0
    assert 0.02 < par_yield(Pu, 10) < 0.06


def test_higher_coupon_is_premium():
    P = _flat(0.04)
    y = par_yield(P, 5)
    assert par_bond_price(P, y + 0.01, 5) > 100.0


def test_validation():
    P = _flat(0.04)
    with pytest.raises(ValueError):
        par_yield(P, 0)
    with pytest.raises(ValueError):
        par_yield(P, 5, freq=0)
