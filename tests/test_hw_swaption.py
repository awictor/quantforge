"""Hull-White swaption (Jamshidian decomposition)."""

import math

import pytest

from quantforge import hw_swaption as sw


def _P0():
    return lambda T: math.exp(-0.03 * T)


EXPIRY = 1.0
PAYS = [2.0, 3.0, 4.0, 5.0]


def _annuity(P0, K):
    prev = EXPIRY
    a = 0.0
    for pt in PAYS:
        a += (pt - prev) * P0(pt)
        prev = pt
    return a


def test_payer_receiver_positive():
    P0 = _P0()
    assert sw(P0, 0.03, 0.1, 0.01, EXPIRY, PAYS, 0.03, True) > 0
    assert sw(P0, 0.03, 0.1, 0.01, EXPIRY, PAYS, 0.03, False) > 0


def test_payer_minus_receiver_is_forward_swap():
    P0 = _P0()
    K = 0.03
    pay = sw(P0, 0.03, 0.1, 0.01, EXPIRY, PAYS, K, True)
    rec = sw(P0, 0.03, 0.1, 0.01, EXPIRY, PAYS, K, False)
    swap = P0(EXPIRY) - P0(PAYS[-1]) - K * _annuity(P0, K)
    assert abs((pay - rec) - swap) < 1e-6


def test_atm_payer_equals_receiver():
    P0 = _P0()
    k_atm = (P0(EXPIRY) - P0(PAYS[-1])) / _annuity(P0, 0.0)
    pa = sw(P0, 0.03, 0.1, 0.01, EXPIRY, PAYS, k_atm, True)
    ra = sw(P0, 0.03, 0.1, 0.01, EXPIRY, PAYS, k_atm, False)
    assert abs(pa - ra) < 1e-6


def test_higher_vol_raises_payer():
    P0 = _P0()
    base = sw(P0, 0.03, 0.1, 0.01, EXPIRY, PAYS, 0.03, True)
    assert sw(P0, 0.03, 0.1, 0.03, EXPIRY, PAYS, 0.03, True) > base


def test_validation():
    P0 = _P0()
    with pytest.raises(ValueError):
        sw(P0, 0.03, 0.1, 0.01, 2.0, [1.0], 0.03)     # pay <= expiry
    with pytest.raises(ValueError):
        sw(P0, 0.03, 0.1, -0.01, EXPIRY, PAYS, 0.03)  # negative sigma
