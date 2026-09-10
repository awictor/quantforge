"""Tests for key-rate (bucketed) DV01."""

import math

import pytest

from quantforge import key_rate_dv01


def test_buckets_sum_to_parallel():
    # A simple bond-like PV: sum of discounted cashflows on a zero curve.
    cashflows = {1.0: 5.0, 2.0: 5.0, 3.0: 105.0}

    def price(curve):
        return sum(cf * math.exp(-curve[t] * t) for t, cf in cashflows.items())

    base = {1.0: 0.02, 2.0: 0.025, 3.0: 0.03}
    kr = key_rate_dv01(price, base)
    assert kr.total_bucketed == pytest.approx(kr.parallel, abs=1e-8)


def test_dv01_negative_for_long_bond():
    # Rates up -> a long bond's PV falls -> negative DV01 in every bucket.
    cashflows = {1.0: 3.0, 2.0: 103.0}

    def price(curve):
        return sum(cf * math.exp(-curve[t] * t) for t, cf in cashflows.items())

    base = {1.0: 0.02, 2.0: 0.03}
    kr = key_rate_dv01(price, base)
    for tenor, dv in kr.buckets.items():
        assert dv < 0


def test_isolated_tenor_only_moves_its_bucket():
    # A cashflow only at t=2 has key-rate DV01 only in the 2y bucket.
    def price(curve):
        return 100.0 * math.exp(-curve[2.0] * 2.0)

    base = {1.0: 0.02, 2.0: 0.03, 3.0: 0.035}
    kr = key_rate_dv01(price, base)
    assert kr.buckets[1.0] == pytest.approx(0.0, abs=1e-10)
    assert kr.buckets[3.0] == pytest.approx(0.0, abs=1e-10)
    assert kr.buckets[2.0] < 0


def test_dv01_magnitude_matches_analytic():
    # Single 3y zero: dPV/dr = -T * PV; DV01 per 1bp = -T * PV * 1e-4.
    def price(curve):
        return 100.0 * math.exp(-curve[3.0] * 3.0)

    base = {3.0: 0.03}
    kr = key_rate_dv01(price, base)
    pv = 100.0 * math.exp(-0.03 * 3.0)
    expected = -3.0 * pv * 1e-4
    assert kr.buckets[3.0] == pytest.approx(expected, rel=1e-4)


def test_scale_normalizes_to_one_bp():
    # A larger bump reports the same 1bp-normalized DV01.
    def price(curve):
        return 100.0 * math.exp(-curve[2.0] * 2.0)

    base = {2.0: 0.03}
    small = key_rate_dv01(price, base, bump=1e-4)
    big = key_rate_dv01(price, base, bump=1e-2)
    assert small.buckets[2.0] == pytest.approx(big.buckets[2.0], rel=1e-3)


def test_one_sided_matches_central_closely():
    cashflows = {1.0: 5.0, 2.0: 105.0}

    def price(curve):
        return sum(cf * math.exp(-curve[t] * t) for t, cf in cashflows.items())

    base = {1.0: 0.02, 2.0: 0.03}
    central = key_rate_dv01(price, base, one_sided=False)
    one = key_rate_dv01(price, base, one_sided=True)
    assert one.parallel == pytest.approx(central.parallel, rel=1e-3)


def test_rejects_bad_bump():
    with pytest.raises(ValueError):
        key_rate_dv01(lambda c: 1.0, {1.0: 0.02}, bump=0.0)
