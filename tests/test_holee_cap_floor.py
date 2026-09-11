"""Ho-Lee cap/floor via the bond-option identity, with swap parity."""

import pytest

from quantforge import holee_cap, holee_floor, holee_caplet
from quantforge.holee import holee_zero_coupon_bond as hz


R0 = 0.03
DATES = [0.5 * i for i in range(1, 9)]
K = 0.03
THETA, SIG = 0.005, 0.01


def test_cap_floor_parity():
    c = holee_cap(R0, DATES, K, THETA, SIG)
    f = holee_floor(R0, DATES, K, THETA, SIG)
    P0, PN = hz(R0, DATES[0], THETA, SIG), hz(R0, DATES[-1], THETA, SIG)
    fixed = K * sum((DATES[i] - DATES[i - 1]) * hz(R0, DATES[i], THETA, SIG)
                    for i in range(1, len(DATES)))
    swap = (P0 - PN) - fixed
    assert (c - f) == pytest.approx(swap, abs=1e-8)


def test_cap_is_sum_of_caplets():
    c = holee_cap(R0, DATES, K, THETA, SIG)
    summed = sum(holee_caplet(R0, DATES[i - 1], DATES[i], K, THETA, SIG)
                 for i in range(1, len(DATES)))
    assert c == pytest.approx(summed, abs=1e-12)


def test_positive_and_monotone():
    assert holee_cap(R0, DATES, K, THETA, SIG) > 0.0
    assert holee_floor(R0, DATES, K, THETA, SIG) > 0.0
    assert (holee_cap(R0, DATES, 0.04, THETA, SIG)
            < holee_cap(R0, DATES, 0.02, THETA, SIG))
    assert (holee_floor(R0, DATES, 0.04, THETA, SIG)
            > holee_floor(R0, DATES, 0.02, THETA, SIG))


def test_too_few_dates_raise():
    with pytest.raises(ValueError):
        holee_cap(R0, [0.5], K, THETA, SIG)
