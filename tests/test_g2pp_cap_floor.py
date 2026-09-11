"""G2++ cap, floor, and caplet/floorlet parity."""

import math

import pytest

from quantforge import g2pp_cap, g2pp_floor, g2pp_caplet, g2pp_floorlet


A, B, SIG, ETA, RHO = 0.1, 0.05, 0.01, 0.008, -0.7
R = 0.03
DS = [(i * 0.5, math.exp(-R * i * 0.5)) for i in range(1, 9)]  # semiannual to 4y
K = 0.03


def test_cap_minus_floor_equals_swap_value():
    c = g2pp_cap(DS, A, B, SIG, ETA, RHO, K)
    f = g2pp_floor(DS, A, B, SIG, ETA, RHO, K)
    P0, PN = DS[0][1], DS[-1][1]
    fixed = K * sum((DS[i][0] - DS[i - 1][0]) * DS[i][1]
                    for i in range(1, len(DS)))
    swap = (P0 - PN) - fixed
    assert (c - f) == pytest.approx(swap, abs=1e-8)


def test_cap_is_sum_of_caplets():
    c = g2pp_cap(DS, A, B, SIG, ETA, RHO, K)
    summed = sum(g2pp_caplet(DS[i - 1][1], DS[i][1], A, B, SIG, ETA, RHO,
                             DS[i - 1][0], DS[i][0], K)
                 for i in range(1, len(DS)))
    assert c == pytest.approx(summed, abs=1e-12)


def test_floor_is_sum_of_floorlets():
    f = g2pp_floor(DS, A, B, SIG, ETA, RHO, K)
    summed = sum(g2pp_floorlet(DS[i - 1][1], DS[i][1], A, B, SIG, ETA, RHO,
                               DS[i - 1][0], DS[i][0], K)
                 for i in range(1, len(DS)))
    assert f == pytest.approx(summed, abs=1e-12)


def test_cap_and_floor_positive():
    assert g2pp_cap(DS, A, B, SIG, ETA, RHO, K) > 0.0
    assert g2pp_floor(DS, A, B, SIG, ETA, RHO, K) > 0.0


def test_higher_strike_lowers_cap_raises_floor():
    lo, hi = 0.02, 0.04
    assert g2pp_cap(DS, A, B, SIG, ETA, RHO, hi) < g2pp_cap(DS, A, B, SIG, ETA,
                                                            RHO, lo)
    assert g2pp_floor(DS, A, B, SIG, ETA, RHO, hi) > g2pp_floor(DS, A, B, SIG,
                                                                ETA, RHO, lo)


def test_too_few_dates_raise():
    with pytest.raises(ValueError):
        g2pp_cap([(0.5, 0.98)], A, B, SIG, ETA, RHO, K)
