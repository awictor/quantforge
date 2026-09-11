"""CDS accrued-premium-on-default (credit.cds_accrual_on_default)."""

import pytest

from quantforge import (
    SurvivalCurve, cds_accrual_on_default, cds_par_spread, risky_annuity,
)
from quantforge.credit import cds_premium_leg


def _curve(h=0.02):
    return SurvivalCurve([1, 3, 5, 10], [h] * 4)


PAY = [0.5 * i for i in range(1, 11)]
R, REC = 0.03, 0.4


def test_accrual_positive_and_small():
    c = _curve()
    a = cds_accrual_on_default(c, PAY, R)
    ann = risky_annuity(c, PAY, R)
    assert a > 0.0
    assert a < 0.02 * ann      # a small fraction of the coupon annuity


def test_accrual_lowers_par_spread():
    c = _curve()
    with_acc = cds_par_spread(c, PAY, R, REC, accrual_on_default=True)
    without = cds_par_spread(c, PAY, R, REC, accrual_on_default=False)
    assert with_acc < without


def test_accrual_raises_premium_leg():
    c = _curve()
    s = cds_par_spread(c, PAY, R, REC)
    assert cds_premium_leg(c, s, PAY, R, accrual_on_default=True) > \
        cds_premium_leg(c, s, PAY, R)


def test_accrual_bounded_by_half_period_default_prob():
    c = _curve()
    a = cds_accrual_on_default(c, PAY, R)
    assert a < 0.5 * (1.0 - c.survival(PAY[-1]))


def test_higher_hazard_more_accrual():
    lo = cds_accrual_on_default(_curve(0.01), PAY, R)
    hi = cds_accrual_on_default(_curve(0.05), PAY, R)
    assert hi > lo
