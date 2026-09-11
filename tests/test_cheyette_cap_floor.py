"""Cheyette cap, floor, and caplet/floorlet parity."""

import math

import pytest

from quantforge import (
    cheyette_cap, cheyette_floor, cheyette_caplet, cheyette_floorlet,
)


KAPPA, SIG = 0.1, 0.01
R = 0.03
DS = [(i * 0.5, math.exp(-R * i * 0.5)) for i in range(1, 9)]
K = 0.03


def test_cap_minus_floor_equals_swap_value():
    c = cheyette_cap(DS, KAPPA, SIG, K)
    f = cheyette_floor(DS, KAPPA, SIG, K)
    P0, PN = DS[0][1], DS[-1][1]
    fixed = K * sum((DS[i][0] - DS[i - 1][0]) * DS[i][1]
                    for i in range(1, len(DS)))
    swap = (P0 - PN) - fixed
    assert (c - f) == pytest.approx(swap, abs=1e-8)


def test_cap_is_sum_of_caplets():
    c = cheyette_cap(DS, KAPPA, SIG, K)
    summed = sum(cheyette_caplet(DS[i - 1][1], DS[i][1], KAPPA, SIG,
                                 DS[i - 1][0], DS[i][0], K)
                 for i in range(1, len(DS)))
    assert c == pytest.approx(summed, abs=1e-12)


def test_floor_is_sum_of_floorlets():
    f = cheyette_floor(DS, KAPPA, SIG, K)
    summed = sum(cheyette_floorlet(DS[i - 1][1], DS[i][1], KAPPA, SIG,
                                   DS[i - 1][0], DS[i][0], K)
                 for i in range(1, len(DS)))
    assert f == pytest.approx(summed, abs=1e-12)


def test_positive_and_monotone():
    assert cheyette_cap(DS, KAPPA, SIG, K) > 0.0
    assert cheyette_floor(DS, KAPPA, SIG, K) > 0.0
    assert cheyette_cap(DS, KAPPA, SIG, 0.04) < cheyette_cap(DS, KAPPA, SIG, 0.02)
    assert cheyette_floor(DS, KAPPA, SIG, 0.04) > cheyette_floor(DS, KAPPA, SIG, 0.02)


def test_too_few_dates_raise():
    with pytest.raises(ValueError):
        cheyette_cap([(0.5, 0.98)], KAPPA, SIG, K)
