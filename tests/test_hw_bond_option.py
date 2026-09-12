"""Hull-White analytic bond option."""

import math

import pytest

from quantforge import hw_bond_option as opt


def _P0():
    return lambda T: math.exp(-0.03 * T)


def test_atm_call_equals_put_and_positive():
    P0 = _P0()
    To, Tb = 1.0, 5.0
    fwd = P0(Tb) / P0(To)
    c = opt(P0, 0.1, 0.01, To, Tb, fwd, True)
    p = opt(P0, 0.1, 0.01, To, Tb, fwd, False)
    assert c > 0 and p > 0
    assert abs(c - p) < 1e-6


def test_put_call_parity():
    P0 = _P0()
    To, Tb, K = 1.0, 5.0, 0.85
    c = opt(P0, 0.1, 0.01, To, Tb, K, True)
    p = opt(P0, 0.1, 0.01, To, Tb, K, False)
    assert abs((c - p) - (P0(Tb) - K * P0(To))) < 1e-9


def test_higher_vol_raises_price():
    P0 = _P0()
    To, Tb = 1.0, 5.0
    fwd = P0(Tb) / P0(To)
    lo = opt(P0, 0.1, 0.01, To, Tb, fwd, True)
    hi = opt(P0, 0.1, 0.03, To, Tb, fwd, True)
    assert hi > lo


def test_deep_itm_approaches_intrinsic():
    P0 = _P0()
    To, Tb, K = 1.0, 5.0, 0.5
    fwd = P0(Tb) / P0(To)
    c = opt(P0, 0.1, 0.01, To, Tb, K, True)
    assert abs(c - (fwd - K) * P0(To)) < 1e-4


def test_zero_vol_atm_worthless():
    P0 = _P0()
    To, Tb = 1.0, 5.0
    fwd = P0(Tb) / P0(To)
    assert opt(P0, 0.1, 0.0, To, Tb, fwd, True) < 1e-9


def test_validation():
    P0 = _P0()
    with pytest.raises(ValueError):
        opt(P0, 0.1, 0.01, 5.0, 1.0, 0.9)        # t_option > t_bond
    with pytest.raises(ValueError):
        opt(P0, 0.1, -0.01, 1.0, 5.0, 0.9)       # negative sigma
