"""Bootstrap a survival curve from par CDS quotes (credit.bootstrap_survival_curve)."""

import pytest

from quantforge import bootstrap_survival_curve, cds_par_spread, SurvivalCurve


MATS = [1, 3, 5, 7, 10]
SPREADS = [0.008, 0.010, 0.012, 0.013, 0.014]
R, REC = 0.03, 0.4


def _repriced(curve, T):
    freq = 4
    pay = [k / freq for k in range(1, T * freq + 1)]
    return cds_par_spread(curve, pay, R, REC, max(1, T * 100))


def test_reprices_all_quotes():
    c = bootstrap_survival_curve(MATS, SPREADS, R, REC)
    for T, s in zip(MATS, SPREADS):
        assert _repriced(c, T) == pytest.approx(s, abs=1e-6)


def test_returns_survival_curve_with_right_pillars():
    c = bootstrap_survival_curve(MATS, SPREADS, R, REC)
    assert isinstance(c, SurvivalCurve)
    assert c.times == MATS
    assert len(c.hazards) == len(MATS)


def test_survival_monotone_decreasing():
    c = bootstrap_survival_curve(MATS, SPREADS, R, REC)
    assert c.survival(1.0) > c.survival(5.0) > c.survival(10.0)


def test_upward_spreads_give_rising_hazards():
    c = bootstrap_survival_curve(MATS, SPREADS, R, REC)
    assert c.hazards[0] < c.hazards[-1]


def test_flat_spreads_give_flat_hazard():
    flat = [0.01] * 5
    c = bootstrap_survival_curve(MATS, flat, R, REC)
    # Nearly constant hazard across pillars for a flat spread term structure.
    assert max(c.hazards) - min(c.hazards) < 0.002


def test_validation():
    with pytest.raises(ValueError):
        bootstrap_survival_curve([1, 3], [0.01], R, REC)
    with pytest.raises(ValueError):
        bootstrap_survival_curve([], [], R, REC)
