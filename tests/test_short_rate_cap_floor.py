"""Vasicek and CIR cap/floor via the bond-option identity, with swap parity."""

import pytest

from quantforge import (
    vasicek_cap, vasicek_floor, vasicek_caplet,
    cir_cap, cir_floor, cir_caplet,
)
from quantforge.vasicek import zero_coupon_bond as vz
from quantforge.cir import cir_zero_coupon_bond as cz


R0 = 0.04
DATES = [0.5 * i for i in range(1, 9)]
K = 0.03
VK, VT, VS = 0.15, 0.05, 0.01
CK, CT, CS = 0.3, 0.04, 0.05


def _swap(df):
    P0, PN = df(DATES[0]), df(DATES[-1])
    fixed = K * sum((DATES[i] - DATES[i - 1]) * df(DATES[i])
                    for i in range(1, len(DATES)))
    return (P0 - PN) - fixed


def test_vasicek_cap_floor_parity():
    c = vasicek_cap(R0, DATES, K, VK, VT, VS)
    f = vasicek_floor(R0, DATES, K, VK, VT, VS)
    swap = _swap(lambda t: vz(R0, t, VK, VT, VS))
    assert (c - f) == pytest.approx(swap, abs=1e-8)


def test_cir_cap_floor_parity():
    c = cir_cap(R0, DATES, K, CK, CT, CS)
    f = cir_floor(R0, DATES, K, CK, CT, CS)
    swap = _swap(lambda t: cz(R0, t, CK, CT, CS))
    assert (c - f) == pytest.approx(swap, abs=1e-8)


def test_vasicek_cap_is_sum_of_caplets():
    c = vasicek_cap(R0, DATES, K, VK, VT, VS)
    summed = sum(vasicek_caplet(R0, DATES[i - 1], DATES[i], K, VK, VT, VS)
                 for i in range(1, len(DATES)))
    assert c == pytest.approx(summed, abs=1e-12)


def test_cir_cap_is_sum_of_caplets():
    c = cir_cap(R0, DATES, K, CK, CT, CS)
    summed = sum(cir_caplet(R0, DATES[i - 1], DATES[i], K, CK, CT, CS)
                 for i in range(1, len(DATES)))
    assert c == pytest.approx(summed, abs=1e-12)


def test_positive_and_monotone():
    assert vasicek_cap(R0, DATES, K, VK, VT, VS) > 0.0
    assert cir_floor(R0, DATES, K, CK, CT, CS) > 0.0
    assert (vasicek_cap(R0, DATES, 0.04, VK, VT, VS)
            < vasicek_cap(R0, DATES, 0.02, VK, VT, VS))


def test_too_few_dates_raise():
    with pytest.raises(ValueError):
        cir_cap(R0, [0.5], K, CK, CT, CS)
