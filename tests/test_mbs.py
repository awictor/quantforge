"""Mortgage-backed security cashflows and prepayment conventions."""

import pytest

from quantforge import (
    monthly_payment, cpr_to_smm, smm_to_cpr, psa_cpr, amortization_schedule,
    mbs_cashflows, weighted_average_life,
    mbs_cashflows_psa, mbs_price, mbs_yield,
    mbs_price_with_spread, mbs_zspread, mbs_effective_duration,
    mbs_effective_convexity,
)


def test_monthly_payment_value():
    assert monthly_payment(300000, 0.05, 360) == pytest.approx(1610.46, abs=0.01)


def test_zero_rate_is_straight_line():
    assert monthly_payment(120000, 0.0, 360) == pytest.approx(120000 / 360)


def test_amortization_pays_off():
    sched = amortization_schedule(300000, 0.05, 360)
    assert sched[-1][3] == pytest.approx(0.0, abs=1e-6)
    assert sum(r[2] for r in sched) == pytest.approx(300000, abs=1e-3)


def test_cpr_smm_round_trip():
    assert smm_to_cpr(cpr_to_smm(0.06)) == pytest.approx(0.06, abs=1e-12)


def test_psa_ramp():
    assert psa_cpr(1, 100) == pytest.approx(0.002)
    assert psa_cpr(30, 100) == pytest.approx(0.06)
    assert psa_cpr(60, 100) == pytest.approx(0.06)  # flat after month 30
    assert psa_cpr(30, 200) == pytest.approx(0.12)


def test_mbs_no_prepay_matches_amortization():
    rows = mbs_cashflows(300000, 0.05, 360, 0.0)
    assert sum(r[4] for r in rows) == pytest.approx(300000, abs=1e-3)


def test_mbs_prepay_principal_sums_to_balance():
    smm = cpr_to_smm(0.06)
    rows = mbs_cashflows(300000, 0.05, 360, smm)
    assert sum(r[4] for r in rows) == pytest.approx(300000, abs=1e-2)


def test_wal_shortens_with_prepayment():
    base = mbs_cashflows(300000, 0.05, 360, 0.0)
    fast = mbs_cashflows(300000, 0.05, 360, cpr_to_smm(0.06))
    assert weighted_average_life(fast, 300000) < weighted_average_life(base, 300000)


def test_zspread_flat_curve_matches_price():
    cf = mbs_cashflows_psa(300000, 0.05, 360, 100)
    flat = [0.05] * len(cf)
    assert mbs_price_with_spread(cf, flat, 0.0) == pytest.approx(
        mbs_price(cf, 0.05), abs=1e-6)


def test_zspread_inverts():
    cf = mbs_cashflows_psa(300000, 0.05, 360, 100)
    flat = [0.05] * len(cf)
    target = mbs_price_with_spread(cf, flat, 0.01)
    assert mbs_zspread(cf, flat, target) == pytest.approx(0.01, abs=1e-8)


def test_zspread_sloped_curve():
    cf = mbs_cashflows_psa(300000, 0.05, 360, 100)
    n = len(cf)
    slope = [0.03 + 0.02 * (i / n) for i in range(n)]
    target = mbs_price_with_spread(cf, slope, 0.005)
    assert mbs_zspread(cf, slope, target) == pytest.approx(0.005, abs=1e-8)


def test_effective_duration_positive_matches_fd():
    cf = mbs_cashflows_psa(300000, 0.05, 360, 100)
    d = mbs_effective_duration(cf, 0.05)
    assert d > 0
    h = 1e-6
    fd = -(mbs_price(cf, 0.05 + h) - mbs_price(cf, 0.05 - h)) / (
        2 * h) / mbs_price(cf, 0.05)
    assert d == pytest.approx(fd, abs=1e-2)


def test_effective_convexity_order():
    cf = mbs_cashflows_psa(300000, 0.05, 360, 100)
    # Convexity is roughly duration^2 in magnitude for a bond-like profile.
    d = mbs_effective_duration(cf, 0.05)
    c = mbs_effective_convexity(cf, 0.05)
    assert 0 < c < 10 * d * d


def test_spread_duration_validation():
    cf = mbs_cashflows_psa(300000, 0.05, 360, 100)
    with pytest.raises(ValueError):
        mbs_zspread(cf, [0.05] * len(cf), -5)
    with pytest.raises(ValueError):
        mbs_price_with_spread(cf, [0.05], 0.0)


def test_psa_zero_reduces_to_no_prepay():
    z = mbs_cashflows_psa(300000, 0.05, 360, 0)
    base = mbs_cashflows(300000, 0.05, 360, 0.0)
    assert sum(r[4] for r in z) == pytest.approx(sum(r[4] for r in base), abs=1e-6)


def test_psa_principal_sums_to_balance():
    rows = mbs_cashflows_psa(300000, 0.05, 360, 100)
    assert sum(r[4] for r in rows) == pytest.approx(300000, abs=1e-2)


def test_faster_psa_shortens_wal():
    w100 = weighted_average_life(mbs_cashflows_psa(300000, 0.05, 360, 100), 300000)
    w300 = weighted_average_life(mbs_cashflows_psa(300000, 0.05, 360, 300), 300000)
    assert w300 < w100


def test_price_monotone_decreasing_in_yield():
    rows = mbs_cashflows_psa(300000, 0.05, 360, 100)
    assert mbs_price(rows, 0.07) < mbs_price(rows, 0.05)


def test_price_at_coupon_is_par():
    rows = mbs_cashflows_psa(300000, 0.05, 360, 100)
    assert mbs_price(rows, 0.05) == pytest.approx(300000, abs=1.0)


def test_yield_inverts_price():
    rows = mbs_cashflows_psa(300000, 0.05, 360, 100)
    p = mbs_price(rows, 0.05)
    assert mbs_yield(rows, p) == pytest.approx(0.05, abs=1e-8)


def test_price_yield_validation():
    with pytest.raises(ValueError):
        mbs_yield(mbs_cashflows_psa(300000, 0.05, 360, 100), -5)


def test_validation():
    with pytest.raises(ValueError):
        cpr_to_smm(1.5)
    with pytest.raises(ValueError):
        monthly_payment(-1, 0.05, 360)
    with pytest.raises(ValueError):
        psa_cpr(0)
