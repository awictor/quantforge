"""Mortgage-backed security cashflows and prepayment conventions."""

import pytest

from quantforge import (
    monthly_payment, cpr_to_smm, smm_to_cpr, psa_cpr, amortization_schedule,
    mbs_cashflows, weighted_average_life,
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


def test_validation():
    with pytest.raises(ValueError):
        cpr_to_smm(1.5)
    with pytest.raises(ValueError):
        monthly_payment(-1, 0.05, 360)
    with pytest.raises(ValueError):
        psa_cpr(0)
