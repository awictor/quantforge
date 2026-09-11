"""CDS risk sensitivities (credit.cds_greeks)."""

import pytest

from quantforge import SurvivalCurve, cds_greeks, cds_par_spread, cds_value


def _curve(h=0.02):
    return SurvivalCurve([1, 3, 5, 10], [h] * 4)


PAY = [0.5 * i for i in range(1, 11)]
R, REC = 0.03, 0.4


def test_value_matches_cds_value():
    c = _curve()
    ps = cds_par_spread(c, PAY, R, REC)
    g = cds_greeks(c, ps, PAY, R, REC)
    assert g["value"] == pytest.approx(cds_value(c, ps, PAY, R, REC), abs=1e-12)


def test_buyer_credit01_positive_recovery01_negative():
    c = _curve()
    ps = cds_par_spread(c, PAY, R, REC)
    g = cds_greeks(c, ps, PAY, R, REC, protection_buyer=True)
    assert g["credit01"] > 0.0     # widening hazards help the buyer
    assert g["recovery01"] < 0.0   # higher recovery hurts protection value


def test_seller_signs_opposite():
    c = _curve()
    ps = cds_par_spread(c, PAY, R, REC)
    buyer = cds_greeks(c, ps, PAY, R, REC, protection_buyer=True)
    seller = cds_greeks(c, ps, PAY, R, REC, protection_buyer=False)
    assert seller["credit01"] == pytest.approx(-buyer["credit01"], abs=1e-12)
    assert seller["value"] == pytest.approx(-buyer["value"], abs=1e-12)


def test_risky_annuity_positive():
    g = cds_greeks(_curve(), 0.012, PAY, R, REC)
    assert g["risky_annuity"] > 0.0


def test_credit01_grows_with_bump():
    c = _curve()
    ps = cds_par_spread(c, PAY, R, REC)
    small = cds_greeks(c, ps, PAY, R, REC, bump=1e-4)["credit01"]
    big = cds_greeks(c, ps, PAY, R, REC, bump=1e-3)["credit01"]
    # Roughly linear: a 10x bump gives ~10x the value change.
    assert big == pytest.approx(10 * small, rel=0.05)
