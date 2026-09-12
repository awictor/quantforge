"""Convertible bond on an equity binomial lattice."""

import pytest

from quantforge import convertible_bond_lattice
from quantforge.equity_comp import conversion_value, straight_bond_floor


S, SIG, FACE, RATIO, CPN, MAT, R = 50.0, 0.3, 1000.0, 20, 0.04, 5, 0.05


def test_above_parity_and_floor():
    v = convertible_bond_lattice(S, SIG, FACE, RATIO, CPN, MAT, R, 200)
    assert v >= conversion_value(S, RATIO) - 1e-6
    assert v >= straight_bond_floor(FACE, CPN, MAT, R, 0.0) - 1e-6


def test_deep_itm_approaches_parity():
    v = convertible_bond_lattice(200, SIG, FACE, RATIO, CPN, MAT, R, 200)
    conv = conversion_value(200, RATIO)
    assert abs(v - conv) / conv < 0.1


def test_callable_below_non_callable():
    v = convertible_bond_lattice(S, SIG, FACE, RATIO, CPN, MAT, R, 200)
    call = convertible_bond_lattice(S, SIG, FACE, RATIO, CPN, MAT, R, 200,
                                    call_price=1100)
    assert call <= v + 1e-6


def test_higher_vol_raises_value():
    v = convertible_bond_lattice(S, SIG, FACE, RATIO, CPN, MAT, R, 200)
    assert convertible_bond_lattice(S, 0.5, FACE, RATIO, CPN, MAT, R, 200) > v


def test_validation():
    with pytest.raises(ValueError):
        convertible_bond_lattice(-1, SIG, FACE, RATIO, CPN, MAT, R)
    with pytest.raises(ValueError):
        convertible_bond_lattice(S, SIG, FACE, RATIO, CPN, 0, R)
