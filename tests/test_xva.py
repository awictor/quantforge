"""Counterparty valuation adjustments (CVA / DVA / BCVA)."""

import pytest

from quantforge import (
    SurvivalCurve, marginal_default_probs, cva, dva, bcva,
)


CURVE = SurvivalCurve([1, 3, 5], [0.02, 0.03, 0.04])
GRID = [1, 2, 3, 4, 5]
EPE = [10, 9, 8, 6, 4]


def test_marginal_default_probs_sum_and_range():
    dq = marginal_default_probs(CURVE, GRID)
    assert sum(dq) == pytest.approx(1 - CURVE.survival(5), abs=1e-12)
    assert all(0.0 <= x <= 1.0 for x in dq)


def test_cva_non_negative():
    assert cva(CURVE, GRID, EPE, 0.03, 0.4) > 0


def test_cva_zero_without_exposure():
    assert cva(CURVE, GRID, [0] * 5, 0.03, 0.4) == pytest.approx(0.0, abs=1e-15)


def test_cva_monotone_in_hazard():
    riskier = SurvivalCurve([1, 3, 5], [0.04, 0.06, 0.08])
    assert cva(riskier, GRID, EPE, 0.03, 0.4) > cva(CURVE, GRID, EPE, 0.03, 0.4)


def test_cva_monotone_in_exposure():
    base = cva(CURVE, GRID, EPE, 0.03, 0.4)
    assert cva(CURVE, GRID, [x * 2 for x in EPE], 0.03, 0.4) > base


def test_cva_lgd_scaling():
    # CVA at zero recovery is exactly twice the CVA at 50% recovery.
    assert cva(CURVE, GRID, EPE, 0.03, 0.0) == pytest.approx(
        2 * cva(CURVE, GRID, EPE, 0.03, 0.5), abs=1e-12)


def test_higher_recovery_lowers_cva():
    assert cva(CURVE, GRID, EPE, 0.03, 0.6) < cva(CURVE, GRID, EPE, 0.03, 0.4)


def test_bcva_is_cva_minus_dva():
    own = SurvivalCurve([1, 3, 5], [0.01, 0.015, 0.02])
    ene = [3, 3, 2, 2, 1]
    expect = cva(CURVE, GRID, EPE, 0.03) - dva(own, GRID, ene, 0.03)
    assert bcva(CURVE, own, GRID, EPE, ene, 0.03) == pytest.approx(expect, abs=1e-15)


def test_bcva_symmetric_is_zero():
    # Identical curves and exposures: CVA and DVA cancel.
    assert bcva(CURVE, CURVE, GRID, EPE, EPE, 0.03) == pytest.approx(0.0, abs=1e-15)


def test_validation():
    with pytest.raises(ValueError):
        cva(CURVE, GRID, EPE, 0.03, 1.5)
    with pytest.raises(ValueError):
        cva(CURVE, [1, 2], EPE, 0.03)
    with pytest.raises(ValueError):
        marginal_default_probs(CURVE, [2, 1])
    with pytest.raises(ValueError):
        cva(CURVE, GRID, [-1, 0, 0, 0, 0], 0.03)
