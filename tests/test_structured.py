"""Structured notes: principal-protected notes and reverse convertibles."""

import math

import pytest

from quantforge import (
    note_zero_coupon_bond, principal_protected_note, reverse_convertible,
    note_embedded_option_value, capped_principal_protected_note,
    reverse_convertible_fair_coupon, buffered_note,
)
from quantforge.bsm import call_price, put_price


S, K, T, R, SIG = 100.0, 100.0, 3.0, 0.04, 0.25


def test_note_zero_coupon_bond():
    assert note_zero_coupon_bond(1000, R, T) == pytest.approx(1000 * math.exp(-R * T))


def test_ppn_above_floor():
    p = principal_protected_note(S, K, T, R, SIG, 1000)
    assert p >= note_zero_coupon_bond(1000, R, T)


def test_ppn_decomposition():
    p = principal_protected_note(S, K, T, R, SIG, 1000)
    ns = 1000 / S
    expect = note_zero_coupon_bond(1000, R, T) + ns * call_price(S, K, T, R, SIG)
    assert p == pytest.approx(expect, abs=1e-9)


def test_ppn_participation_raises_value():
    base = principal_protected_note(S, K, T, R, SIG, 1000)
    assert principal_protected_note(S, K, T, R, SIG, 1000, participation=2.0) > base


def test_embedded_option_positive_for_ppn():
    p = principal_protected_note(S, K, T, R, SIG, 1000)
    emb = note_embedded_option_value(p, 1000, R, T)
    assert emb > 0
    assert emb == pytest.approx((1000 / S) * call_price(S, K, T, R, SIG), abs=1e-9)


def test_reverse_convertible_decomposition():
    rc = reverse_convertible(S, K, T, R, SIG, 1000, 0.08)
    ns = 1000 / K
    bc = 1000 * (1 + 0.08 * T) * math.exp(-R * T)
    assert rc == pytest.approx(bc - ns * put_price(S, K, T, R, SIG), abs=1e-9)


def test_reverse_convertible_coupon_raises_value():
    base = reverse_convertible(S, K, T, R, SIG, 1000, 0.08)
    assert reverse_convertible(S, K, T, R, SIG, 1000, 0.12) > base


def test_capped_ppn_below_uncapped_above_floor():
    capped = capped_principal_protected_note(S, K, 130, T, R, SIG, 1000)
    uncapped = principal_protected_note(S, K, T, R, SIG, 1000)
    assert capped <= uncapped + 1e-9
    assert capped >= note_zero_coupon_bond(1000, R, T)


def test_higher_cap_approaches_uncapped():
    low = capped_principal_protected_note(S, K, 130, T, R, SIG, 1000)
    high = capped_principal_protected_note(S, K, 200, T, R, SIG, 1000)
    assert high > low


def test_fair_coupon_prices_at_par():
    cpn = reverse_convertible_fair_coupon(S, K, T, R, SIG, 1000)
    assert cpn > 0
    assert reverse_convertible(S, K, T, R, SIG, 1000, cpn) == pytest.approx(1000, abs=1e-6)


def test_buffer_raises_value_monotonically():
    b0 = buffered_note(S, K, 0.0, T, R, SIG, 1000)
    b10 = buffered_note(S, K, 0.10, T, R, SIG, 1000)
    b20 = buffered_note(S, K, 0.20, T, R, SIG, 1000)
    assert b0 < b10 < b20


def test_capped_buffer_validation():
    with pytest.raises(ValueError):
        capped_principal_protected_note(S, K, 90, T, R, SIG, 1000)
    with pytest.raises(ValueError):
        buffered_note(S, K, 1.5, T, R, SIG, 1000)


def test_validation():
    with pytest.raises(ValueError):
        principal_protected_note(S, K, T, R, SIG, 1000, participation=-1)
    with pytest.raises(ValueError):
        reverse_convertible(S, K, T, R, SIG, 1000, -0.01)
    with pytest.raises(ValueError):
        note_zero_coupon_bond(-1, R, T)
