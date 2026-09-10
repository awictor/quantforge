"""Tests for the scenario / stress grid.

We anchor the grid P&L to first-order Greeks (small shocks) and to intuitive
worst-case behavior (a long straddle's worst case is no move; a short
straddle's worst case is a large move).
"""

import pytest

from quantforge import (
    Contract, price_book, stress_grid, spot_ladder, greeks, OptionType,
)


def _call_book(qty=1):
    return [Contract(S=100, K=100, t=0.5, r=0.04, sigma=0.25,
                     option_type="call", qty=qty, multiplier=1)]


def test_zero_shock_is_zero_pnl():
    grid = stress_grid(_call_book(), [0.0], [0.0])
    assert grid.pnl[0][0] == pytest.approx(0.0, abs=1e-12)


def test_small_spot_shock_matches_delta():
    # Absolute +0.01 spot move, vol unchanged: P&L ~ delta * dS.
    contracts = _call_book()
    grid = stress_grid(contracts, [0.01], [0.0], relative=False)
    d = greeks(100, 100, 0.5, 0.04, 0.25, OptionType.CALL).delta
    assert grid.pnl[0][0] == pytest.approx(d * 0.01, abs=1e-4)


def test_small_vol_shock_matches_vega():
    # Absolute +0.001 vol move, spot unchanged: P&L ~ vega * dsigma.
    contracts = _call_book()
    grid = stress_grid(contracts, [0.0], [0.001], relative=False)
    v = greeks(100, 100, 0.5, 0.04, 0.25, OptionType.CALL).vega
    assert grid.pnl[0][0] == pytest.approx(v * 0.001, abs=1e-3)


def test_long_call_gains_when_spot_rises():
    grid = stress_grid(_call_book(), [-0.1, 0.0, 0.1], [0.0], relative=True)
    down, flat, up = grid.pnl[0][0], grid.pnl[1][0], grid.pnl[2][0]
    assert down < flat < up


def test_long_straddle_worst_case_is_no_move():
    straddle = [
        Contract(S=100, K=100, t=0.25, r=0.02, sigma=0.3, option_type="call", qty=1),
        Contract(S=100, K=100, t=0.25, r=0.02, sigma=0.3, option_type="put", qty=1),
    ]
    shocks = [-0.2, -0.1, 0.0, 0.1, 0.2]
    grid = stress_grid(straddle, shocks, [0.0], relative=True)
    ss, vs, worst_pnl = grid.worst_case()
    # A long straddle loses most (time/vol held) when the underlying doesn't move.
    assert ss == pytest.approx(0.0)


def test_short_straddle_worst_case_is_big_move():
    straddle = [
        Contract(S=100, K=100, t=0.25, r=0.02, sigma=0.3, option_type="call", qty=-1),
        Contract(S=100, K=100, t=0.25, r=0.02, sigma=0.3, option_type="put", qty=-1),
    ]
    shocks = [-0.3, -0.15, 0.0, 0.15, 0.3]
    grid = stress_grid(straddle, shocks, [0.0], relative=True)
    ss, vs, worst_pnl = grid.worst_case()
    # Short gamma: worst loss at the largest move on the grid.
    assert abs(ss) == pytest.approx(0.3)
    assert worst_pnl < 0


def test_vega_shock_hurts_short_vol_book():
    # Short a call, then bump vol up: value of the short rises against us -> loss.
    contracts = _call_book(qty=-1)
    grid = stress_grid(contracts, [0.0], [0.5], relative=True)  # +50% vol
    assert grid.pnl[0][0] < 0


def test_grid_dimensions_and_rows():
    grid = stress_grid(_call_book(), [-0.1, 0.0, 0.1], [-0.2, 0.0, 0.2])
    assert len(grid.pnl) == 3
    assert len(grid.pnl[0]) == 3
    rows = list(grid.as_rows())
    assert len(rows) == 3
    assert set(rows[0][1].keys()) == {-0.2, 0.0, 0.2}


def test_best_case_beats_worst_case():
    grid = stress_grid(_call_book(), [-0.2, 0.0, 0.2], [-0.1, 0.0, 0.1])
    _, _, worst = grid.worst_case()
    _, _, best = grid.best_case()
    assert best >= worst


def test_spot_ladder_matches_grid_column():
    contracts = _call_book()
    ladder = spot_ladder(contracts, [-0.1, 0.0, 0.1])
    grid = stress_grid(contracts, [-0.1, 0.0, 0.1], [0.0])
    for (ss, pnl), i in zip(ladder, range(3)):
        assert pnl == pytest.approx(grid.pnl[i][0])
