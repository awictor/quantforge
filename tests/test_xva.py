"""Counterparty valuation adjustments (CVA / DVA / BCVA)."""

import pytest

import math

from quantforge import (
    SurvivalCurve, marginal_default_probs, cva, dva, bcva,
    swap_expected_exposure, fva,
    swap_potential_future_exposure, wrong_way_cva,
)
from quantforge.mathfns import norm_ppf


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


EE_GRID = [0, 1, 2, 3, 4, 5]


def test_swap_exposure_zero_at_ends():
    ee = swap_expected_exposure(1e6, 0.01, 5.0, EE_GRID)
    assert ee[0] == pytest.approx(0.0, abs=1e-12)
    assert ee[-1] == pytest.approx(0.0, abs=1e-12)


def test_swap_exposure_humped_and_nonneg():
    ee = swap_expected_exposure(1e6, 0.01, 5.0, EE_GRID)
    assert all(x >= 0 for x in ee)
    peak = ee.index(max(ee))
    assert 0 < peak < len(ee) - 1


def test_swap_exposure_formula():
    ee = swap_expected_exposure(1e6, 0.01, 5.0, EE_GRID)
    std = 1e6 * 0.01 * math.sqrt(2) * (5 - 2) / 5
    assert ee[2] == pytest.approx(std / math.sqrt(2 * math.pi), abs=1e-6)


def test_fva_proportional_to_spread():
    ee = swap_expected_exposure(1e6, 0.01, 5.0, EE_GRID)
    assert fva(EE_GRID, ee, 0.010, 0.03) == pytest.approx(
        2 * fva(EE_GRID, ee, 0.005, 0.03), abs=1e-9)


def test_fva_zero_at_zero_spread():
    ee = swap_expected_exposure(1e6, 0.01, 5.0, EE_GRID)
    assert fva(EE_GRID, ee, 0.0, 0.03) == pytest.approx(0.0, abs=1e-15)


def test_fva_survival_weighting_reduces():
    ee = swap_expected_exposure(1e6, 0.01, 5.0, EE_GRID)
    curve = SurvivalCurve([5], [0.05])
    assert fva(EE_GRID, ee, 0.005, 0.03, curve.survival) < fva(EE_GRID, ee, 0.005, 0.03)


def test_pfe_above_epe():
    p = swap_potential_future_exposure(1e6, 0.01, 5.0, GRID, 0.95)
    e = swap_expected_exposure(1e6, 0.01, 5.0, GRID)
    assert all(p[i] >= e[i] - 1e-9 for i in range(len(GRID)))


def test_pfe_formula():
    p = swap_potential_future_exposure(1e6, 0.01, 5.0, GRID, 0.95)
    std = 1e6 * 0.01 * math.sqrt(2) * (5 - 2) / 5
    assert p[1] == pytest.approx(std * norm_ppf(0.95), abs=1e-6)


def test_pfe_quantile_monotone():
    p95 = swap_potential_future_exposure(1e6, 0.01, 5.0, GRID, 0.95)
    p99 = swap_potential_future_exposure(1e6, 0.01, 5.0, GRID, 0.99)
    assert all(p99[i] >= p95[i] - 1e-9 for i in range(len(GRID)))


def test_wrong_way_reduces_to_cva_at_zero_alpha():
    assert wrong_way_cva(CURVE, GRID, EPE, 0.03, 0.4, 0.0) == pytest.approx(
        cva(CURVE, GRID, EPE, 0.03, 0.4), abs=1e-12)


def test_wrong_way_raises_cva_for_rising_exposure():
    ee_up = [4, 6, 8, 9, 10]
    base = cva(CURVE, GRID, ee_up, 0.03)
    assert wrong_way_cva(CURVE, GRID, ee_up, 0.03, 0.4, 0.5) > base
    assert wrong_way_cva(CURVE, GRID, ee_up, 0.03, 0.4, -0.5) < base  # right-way


def test_pfe_wrong_way_validation():
    with pytest.raises(ValueError):
        swap_potential_future_exposure(1e6, 0.01, 5.0, GRID, 0.3)


def test_swap_exposure_fva_validation():
    with pytest.raises(ValueError):
        swap_expected_exposure(1e6, 0.01, 5.0, [6])
    with pytest.raises(ValueError):
        swap_expected_exposure(1e6, 0.01, 0, EE_GRID)


def test_validation():
    with pytest.raises(ValueError):
        cva(CURVE, GRID, EPE, 0.03, 1.5)
    with pytest.raises(ValueError):
        cva(CURVE, [1, 2], EPE, 0.03)
    with pytest.raises(ValueError):
        marginal_default_probs(CURVE, [2, 1])
    with pytest.raises(ValueError):
        cva(CURVE, GRID, [-1, 0, 0, 0, 0], 0.03)
