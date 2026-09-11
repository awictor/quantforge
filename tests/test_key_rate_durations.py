"""Key-rate and effective durations against a zero curve (bondmath)."""

import math

import pytest

from quantforge import (
    bond_cashflows, key_rate_durations, effective_duration_from_curve,
    macaulay_duration, price_from_curve, DiscountCurve,
)


CF = bond_cashflows(100.0, 0.05, 5.0, freq=2)
PTS = [0.5, 1.0, 2.0, 3.0, 5.0]
ZR = [0.03, 0.032, 0.035, 0.037, 0.04]


def test_sum_of_krds_equals_effective_duration():
    krds = key_rate_durations(CF, PTS, ZR)
    assert sum(krds) == pytest.approx(effective_duration_from_curve(CF, PTS, ZR),
                                      abs=1e-2)


def test_flat_curve_sum_matches_macaulay():
    flat = [0.04] * 5
    krds = key_rate_durations(CF, PTS, flat)
    assert sum(krds) == pytest.approx(macaulay_duration(CF, 0.04), abs=1e-2)


def test_longest_pillar_carries_most():
    # A bullet bond's cashflows cluster at the final pillar (face + last coupon).
    krds = key_rate_durations(CF, PTS, ZR)
    assert krds[-1] == max(krds)
    assert all(k >= 0 for k in krds)


def test_price_from_curve_matches_manual():
    curve = DiscountCurve.from_zero_rates(PTS, ZR)
    manual = sum(cf * curve.df(t) for t, cf in CF)
    assert price_from_curve(CF, curve) == pytest.approx(manual, abs=1e-12)


def test_price_from_curve_accepts_plain_callable():
    p = price_from_curve(CF, lambda t: math.exp(-0.04 * t))
    assert p > 0.0


def test_validation():
    with pytest.raises(ValueError):
        key_rate_durations(CF, PTS, [0.03])
    with pytest.raises(ValueError):
        effective_duration_from_curve(CF, PTS, [0.03])
