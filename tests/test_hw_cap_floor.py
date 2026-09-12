"""Hull-White cap/floor."""

import math

import pytest

from quantforge import hw_cap, hw_floor, hw_caplet


def _P0():
    return lambda T: math.exp(-0.03 * T)


DATES = [1.0, 2.0, 3.0, 4.0, 5.0]


def test_cap_floor_positive():
    P0 = _P0()
    assert hw_cap(P0, 0.1, 0.01, DATES, 0.03) > 0
    assert hw_floor(P0, 0.1, 0.01, DATES, 0.03) > 0


def test_cap_is_sum_of_caplets():
    P0 = _P0()
    cap = hw_cap(P0, 0.1, 0.01, DATES, 0.03)
    parts = sum(hw_caplet(P0, 0.1, 0.01, DATES[i], DATES[i + 1], 0.03)
                for i in range(len(DATES) - 1))
    assert abs(cap - parts) < 1e-12


def test_cap_minus_floor_is_swap_value():
    P0 = _P0()
    cap = hw_cap(P0, 0.1, 0.01, DATES, 0.03)
    floor = hw_floor(P0, 0.1, 0.01, DATES, 0.03)
    swap = 0.0
    for i in range(len(DATES) - 1):
        tau = DATES[i + 1] - DATES[i]
        fwd = (P0(DATES[i]) / P0(DATES[i + 1]) - 1) / tau
        swap += tau * (fwd - 0.03) * P0(DATES[i + 1])
    assert abs((cap - floor) - swap) < 1e-9


def test_vol_and_strike_monotonic():
    P0 = _P0()
    base = hw_cap(P0, 0.1, 0.01, DATES, 0.03)
    assert hw_cap(P0, 0.1, 0.03, DATES, 0.03) > base    # more vol
    assert hw_cap(P0, 0.1, 0.01, DATES, 0.05) < base    # higher strike


def test_validation():
    P0 = _P0()
    with pytest.raises(ValueError):
        hw_cap(P0, 0.1, 0.01, [1.0], 0.03)              # one date
    with pytest.raises(ValueError):
        hw_caplet(P0, 0.1, 0.01, 2.0, 1.0, 0.03)        # pay < reset
