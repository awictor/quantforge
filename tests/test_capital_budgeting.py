"""Capital-budgeting metrics: NPV, IRR, MIRR, payback, PI."""

import pytest

from quantforge import npv, irr, profitability_index, payback_period, mirr


CF = [-1000, 300, 400, 500, 300]


def test_npv_zero_at_irr():
    r = irr(CF)
    assert npv(r, CF) == pytest.approx(0.0, abs=1e-8)


def test_irr_recovers_known_rate():
    ann = sum(1 / 1.08 ** t for t in range(1, 6))
    F = 1000 / ann
    flows = [-1000] + [F] * 5
    assert irr(flows) == pytest.approx(0.08, abs=1e-6)


def test_profitability_index_at_and_below_irr():
    r = irr(CF)
    assert profitability_index(r, CF) == pytest.approx(1.0, abs=1e-6)
    assert profitability_index(0.05, CF) > 1


def test_payback_period():
    p = payback_period(CF)
    assert 2 < p < 3


def test_mirr_in_bounds():
    m = mirr(CF, 0.05, 0.08)
    assert 0 < m < 1


def test_npv_monotone_in_rate():
    assert npv(0.05, CF) > 0
    assert npv(0.30, CF) < 0


def test_validation():
    with pytest.raises(ValueError):
        irr([-1000, -500])
    with pytest.raises(ValueError):
        profitability_index(0.1, [1000, -300])   # positive outlay
