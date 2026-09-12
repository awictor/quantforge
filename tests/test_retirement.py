"""Retirement drawdown: withdrawals, depletion, glide paths."""

import pytest

from quantforge import (
    portfolio_depletion_years, sustainable_withdrawal, withdrawal_balance_path,
    glide_path_equity_weight, withdrawal_stream_pv, ruin_probability_mc,
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


def test_withdrawal_pv_matches_annuity():
    p = withdrawal_stream_pv(40000, 0.03, 30)
    manual = sum(40000 / (1.03) ** k for k in range(30))
    assert p == pytest.approx(manual, abs=1e-6)


def test_withdrawal_pv_growth_equals_discount():
    assert withdrawal_stream_pv(40000, 0.03, 30, 0.03) == pytest.approx(40000 * 30)


def test_withdrawal_pv_monotonicity():
    base = withdrawal_stream_pv(40000, 0.03, 30)
    assert withdrawal_stream_pv(50000, 0.03, 30) > base
    assert withdrawal_stream_pv(40000, 0.03, 40) > base
    assert withdrawal_stream_pv(40000, 0.05, 30) < base


def test_ruin_rises_with_withdrawal():
    hi = ruin_probability_mc(1e6, 70000, 0.04, 0.12, 30, n_paths=4000)
    lo = ruin_probability_mc(1e6, 30000, 0.04, 0.12, 30, n_paths=4000)
    assert hi > lo


def test_low_withdrawal_near_zero_ruin():
    assert ruin_probability_mc(1e6, 10000, 0.05, 0.10, 30, n_paths=4000) < 0.05


def test_ruin_rises_with_volatility():
    assert ruin_probability_mc(1e6, 50000, 0.04, 0.20, 30, n_paths=4000) > \
        ruin_probability_mc(1e6, 50000, 0.04, 0.05, 30, n_paths=4000)


def test_ruin_deterministic():
    a = ruin_probability_mc(1e6, 50000, 0.04, 0.12, 30, n_paths=2000)
    b = ruin_probability_mc(1e6, 50000, 0.04, 0.12, 30, n_paths=2000)
    assert a == b


def test_pv_ruin_validation():
    with pytest.raises(ValueError):
        withdrawal_stream_pv(-1, 0.03, 30)
    with pytest.raises(ValueError):
        ruin_probability_mc(-1, 40000, 0.04, 0.12, 30)


def test_validation():
    with pytest.raises(ValueError):
        portfolio_depletion_years(-1, 40000, 0.03)
    with pytest.raises(ValueError):
        sustainable_withdrawal(1e6, 0.03, 0)
    with pytest.raises(ValueError):
        glide_path_equity_weight(10, 0, 0.9, 0.3)
