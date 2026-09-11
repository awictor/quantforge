"""Reduced-form credit: survival curve and CDS pricing (credit module)."""

import math

import pytest

from quantforge import (
    SurvivalCurve, risky_annuity, cds_par_spread, cds_value,
)


def _flat(h=0.02):
    return SurvivalCurve([1, 3, 5, 10], [h] * 4)


PAY = [0.5 * i for i in range(1, 11)]  # 5y semiannual
R, REC = 0.03, 0.4


def test_survival_starts_at_one_and_decreases():
    c = _flat()
    assert c.survival(0.0) == 1.0
    assert c.survival(1.0) > c.survival(5.0) > c.survival(10.0)


def test_survival_matches_flat_hazard():
    c = _flat(0.03)
    assert c.survival(4.0) == pytest.approx(math.exp(-0.03 * 4.0), abs=1e-12)


def test_par_spread_near_credit_triangle():
    c = _flat(0.02)
    ps = cds_par_spread(c, PAY, R, REC)
    # h * (1 - recovery), up to discounting/timing effects.
    assert ps == pytest.approx(0.02 * (1 - REC), rel=0.05)


def test_cds_value_zero_at_par():
    c = _flat(0.02)
    ps = cds_par_spread(c, PAY, R, REC)
    assert cds_value(c, ps, PAY, R, REC) == pytest.approx(0.0, abs=1e-6)


def test_protection_buyer_profits_below_par():
    c = _flat(0.02)
    ps = cds_par_spread(c, PAY, R, REC)
    assert cds_value(c, ps * 0.5, PAY, R, REC, protection_buyer=True) > 0.0
    assert cds_value(c, ps * 0.5, PAY, R, REC, protection_buyer=False) < 0.0


def test_higher_hazard_widens_par_spread():
    lo = cds_par_spread(_flat(0.01), PAY, R, REC)
    hi = cds_par_spread(_flat(0.04), PAY, R, REC)
    assert hi > lo


def test_risky_annuity_below_riskfree_annuity():
    c = _flat(0.02)
    risky = risky_annuity(c, PAY, R)
    riskfree = sum(0.5 * math.exp(-R * t) for t in PAY)
    assert risky < riskfree


def test_validation():
    with pytest.raises(ValueError):
        SurvivalCurve([1, 1], [0.02, 0.02])   # non-increasing times
    with pytest.raises(ValueError):
        SurvivalCurve([1, 2], [0.02])          # length mismatch
    with pytest.raises(ValueError):
        SurvivalCurve([1, 2], [0.02, -0.01])   # negative hazard
