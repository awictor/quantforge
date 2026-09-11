"""Overnight compounded / simple-average rates (rates module)."""

import pytest

from quantforge import compounded_overnight_rate, simple_average_rate


def test_compounded_above_simple_for_positive():
    n = 90
    fix = [0.05] * n
    acc = [1 / 360] * n
    assert compounded_overnight_rate(fix, acc) > simple_average_rate(fix, acc)


def test_simple_average_of_flat_is_the_rate():
    assert simple_average_rate([0.05] * 30, [1 / 360] * 30) == pytest.approx(
        0.05, abs=1e-12)


def test_single_fixing_both_equal_rate():
    assert compounded_overnight_rate([0.05], [1 / 360]) == pytest.approx(0.05, abs=1e-9)
    assert simple_average_rate([0.05], [1 / 360]) == pytest.approx(0.05, abs=1e-9)


def test_simple_average_weighted():
    fix = [0.04, 0.06]
    acc = [1 / 360, 3 / 360]
    # weighted average: (0.04*1 + 0.06*3)/4 = 0.055
    assert simple_average_rate(fix, acc) == pytest.approx(0.055, abs=1e-12)


def test_compounded_matches_manual():
    fix = [0.04, 0.05, 0.06]
    acc = [1 / 360] * 3
    growth = 1.0
    for r, tau in zip(fix, acc):
        growth *= (1 + r * tau)
    expected = (growth - 1.0) / sum(acc)
    assert compounded_overnight_rate(fix, acc) == pytest.approx(expected, abs=1e-12)


def test_validation():
    with pytest.raises(ValueError):
        compounded_overnight_rate([0.05], [1 / 360, 1 / 360])
    with pytest.raises(ValueError):
        simple_average_rate([], [])
