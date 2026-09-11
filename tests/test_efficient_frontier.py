"""Efficient frontier and target-return optimization (portopt module)."""

import pytest

from quantforge import (
    target_return_weights, efficient_frontier, min_variance_weights,
    portfolio_variance, portfolio_return,
)


MU = [0.08, 0.05, 0.12]
COV = [[0.04, 0.01, 0.0], [0.01, 0.09, 0.02], [0.0, 0.02, 0.16]]


@pytest.mark.parametrize("tgt", [0.06, 0.09, 0.11])
def test_hits_target_and_sums_to_one(tgt):
    w = target_return_weights(MU, COV, tgt)
    assert portfolio_return(w, MU) == pytest.approx(tgt, abs=1e-9)
    assert sum(w) == pytest.approx(1.0, abs=1e-8)


def test_recovers_min_variance_at_its_own_return():
    wmv = min_variance_weights(COV)
    rmv = portfolio_return(wmv, MU)
    wt = target_return_weights(MU, COV, rmv)
    assert wt == pytest.approx(wmv, abs=1e-8)


def test_frontier_minimum_at_min_variance_return():
    wmv = min_variance_weights(COV)
    rmv = portfolio_return(wmv, MU)
    fr = efficient_frontier(MU, COV, [rmv - 0.03, rmv, rmv + 0.03])
    stds = [s for _, s in fr]
    assert stds[1] == min(stds)   # the min-variance return has the lowest std


def test_frontier_reports_requested_returns():
    fr = efficient_frontier(MU, COV, [0.06, 0.10])
    assert [t for t, _ in fr] == [0.06, 0.10]
    assert all(s > 0 for _, s in fr)


def test_target_variance_above_min_variance():
    wmv = min_variance_weights(COV)
    vmin = portfolio_variance(wmv, COV)
    w = target_return_weights(MU, COV, 0.11)
    assert portfolio_variance(w, COV) >= vmin - 1e-12


def test_validation():
    with pytest.raises(ValueError):
        target_return_weights([0.05], COV, 0.06)   # length mismatch
