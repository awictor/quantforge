"""Portfolio rebalancing: drift, turnover, no-trade bands, cost drag."""

import pytest

from quantforge import (
    drift_weights, turnover, transaction_cost, no_trade_band_rebalance,
)


def test_drift_sums_and_winner_gains():
    d = drift_weights([0.5, 0.5], [0.20, -0.10])
    assert sum(d) == pytest.approx(1.0)
    assert d[0] > 0.5


def test_turnover_zero_on_target():
    assert turnover([0.5, 0.5], [0.5, 0.5]) == 0


def test_turnover_formula_and_full_swap():
    assert turnover([0.6, 0.4], [0.5, 0.5]) == pytest.approx(0.1)
    assert turnover([1.0, 0.0], [0.0, 1.0]) == pytest.approx(1.0)


def test_transaction_cost_formula():
    assert transaction_cost([0.6, 0.4], [0.5, 0.5], 10) == pytest.approx(
        2 * 0.1 * 10 / 1e4)


def test_no_trade_band_suppresses_small_drift():
    r = no_trade_band_rebalance([0.52, 0.48], [0.5, 0.5], 0.05)
    assert r[0] == pytest.approx(0.52)   # within band, untouched (renormalized)


def test_no_trade_band_trades_large_drift():
    r = no_trade_band_rebalance([0.7, 0.3], [0.5, 0.5], 0.05)
    assert r[0] == pytest.approx(0.5)
    assert sum(r) == pytest.approx(1.0)


def test_validation():
    with pytest.raises(ValueError):
        turnover([0.5], [0.5, 0.5])
    with pytest.raises(ValueError):
        transaction_cost([0.5, 0.5], [0.5, 0.5], -1)
