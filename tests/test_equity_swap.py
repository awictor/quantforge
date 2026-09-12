"""Equity swaps: total-return and dividend swaps."""

import pytest

from quantforge import (
    total_return_leg, financing_leg, total_return_swap_value, trs_fair_spread,
    dividend_swap_fair_strike, dividend_swap_value,
)


DIVS = [2.0, 2.1, 2.2]
DFS = [0.97, 0.94, 0.91]


def test_total_return_leg():
    assert total_return_leg(1e6, 100, 110, 3) == pytest.approx(1e6 * (13 / 100))


def test_financing_leg():
    assert financing_leg(1e6, 0.03, 0.005, 1.0) == pytest.approx(1e6 * 0.035)


def test_trs_value_is_equity_minus_financing():
    assert total_return_swap_value(1e6, 100, 110, 3, 0.03, 0.005, 1.0) == \
        pytest.approx(1e6 * 0.13 - 1e6 * 0.035)


def test_fair_spread_zeroes_swap():
    sp = trs_fair_spread(100, 110, 3, 0.03, 1.0)
    assert total_return_swap_value(1e6, 100, 110, 3, 0.03, sp, 1.0) == \
        pytest.approx(0.0, abs=1e-6)


def test_dividend_fair_strike_is_pv():
    K = dividend_swap_fair_strike(DIVS, DFS)
    assert K == pytest.approx(sum(DIVS[i] * DFS[i] for i in range(3)))


def test_dividend_swap_zero_at_fair_strike():
    K = dividend_swap_fair_strike(DIVS, DFS)
    assert dividend_swap_value(DIVS, K, DFS) == pytest.approx(0.0, abs=1e-9)


def test_dividend_swap_gains_when_realized_exceeds():
    K = dividend_swap_fair_strike(DIVS, DFS)
    assert dividend_swap_value([3, 3, 3], K, DFS) > 0


def test_validation():
    with pytest.raises(ValueError):
        total_return_leg(1e6, 0, 110, 3)
    with pytest.raises(ValueError):
        dividend_swap_fair_strike([2], [0.9, 0.8])
    with pytest.raises(ValueError):
        trs_fair_spread(100, 110, 3, 0.03, 0)
