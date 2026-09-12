"""Retirement drawdown: withdrawals, depletion, glide paths."""

import pytest

from quantforge import (
    portfolio_depletion_years, sustainable_withdrawal, withdrawal_balance_path,
    glide_path_equity_weight,
)


def test_depletion_inverts_sustainable_withdrawal():
    W = sustainable_withdrawal(1e6, 0.03, 30)
    assert portfolio_depletion_years(1e6, W, 0.03) == pytest.approx(30, abs=1e-6)


def test_higher_return_extends_horizon():
    assert portfolio_depletion_years(1e6, 50000, 0.05) > \
        portfolio_depletion_years(1e6, 50000, 0.02)


def test_withdrawal_below_interest_is_infinite():
    assert portfolio_depletion_years(1e6, 20000, 0.03) == float("inf")


def test_zero_return_is_balance_over_withdrawal():
    assert portfolio_depletion_years(1e6, 40000, 0.0) == pytest.approx(25)


def test_sustainable_withdrawal_zero_return():
    assert sustainable_withdrawal(1e6, 0.0, 25) == pytest.approx(40000)


def test_inflation_indexed_path_declines():
    path = withdrawal_balance_path(1e6, 40000, 0.05, 0.03, 40)
    assert len(path) == 40
    assert path[-1] < path[0]
    assert min(path) >= 0.0


def test_glide_path_monotone_and_clamped():
    ws = [glide_path_equity_weight(y, 30, 0.9, 0.3) for y in (30, 20, 10, 5, 0)]
    assert all(ws[i] >= ws[i + 1] for i in range(len(ws) - 1))
    assert glide_path_equity_weight(40, 30, 0.9, 0.3) == pytest.approx(0.9)
    assert glide_path_equity_weight(-5, 30, 0.9, 0.3) == pytest.approx(0.3)


def test_validation():
    with pytest.raises(ValueError):
        portfolio_depletion_years(-1, 40000, 0.03)
    with pytest.raises(ValueError):
        sustainable_withdrawal(1e6, 0.03, 0)
    with pytest.raises(ValueError):
        glide_path_equity_weight(10, 0, 0.9, 0.3)
