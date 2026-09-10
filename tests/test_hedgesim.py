"""Tests for the discrete delta-hedge P&L simulator.

The Black-Scholes replication argument says a continuously delta-hedged short
option breaks even in expectation; discrete rehedging leaves a zero-mean error
whose standard deviation shrinks like 1/sqrt(n_steps) (Boyle-Emanuel).
"""

import math

import pytest

from quantforge import simulate_delta_hedge, OptionType


@pytest.mark.slow
def test_mean_pnl_near_zero():
    # Hedging at the true vol: expected P&L ~ 0 (premium replicates payoff).
    res = simulate_delta_hedge(100, 100, 1.0, 0.05, 0.2, OptionType.CALL,
                               n_steps=50, n_paths=30_000, seed=1)
    # Mean should be within a few Monte Carlo standard errors of zero.
    se = res.std_pnl / math.sqrt(res.n_paths)
    assert abs(res.mean_pnl) < 5 * se


@pytest.mark.slow
def test_error_shrinks_with_more_rehedges():
    # Doubling steps ~ divides the hedging-error std by sqrt(2).
    coarse = simulate_delta_hedge(100, 100, 1.0, 0.05, 0.2, OptionType.CALL,
                                  n_steps=25, n_paths=20_000, seed=7)
    fine = simulate_delta_hedge(100, 100, 1.0, 0.05, 0.2, OptionType.CALL,
                                n_steps=100, n_paths=20_000, seed=7)
    # 4x steps -> ~2x tighter. Allow slack for MC noise.
    ratio = coarse.std_pnl / fine.std_pnl
    assert 1.6 < ratio < 2.5


@pytest.mark.slow
def test_put_also_replicates():
    res = simulate_delta_hedge(100, 100, 1.0, 0.05, 0.25, OptionType.PUT,
                               n_steps=50, n_paths=30_000, seed=3)
    se = res.std_pnl / math.sqrt(res.n_paths)
    assert abs(res.mean_pnl) < 5 * se


def test_hedging_at_low_vol_loses_when_realized_is_high():
    # Short an option priced/hedged at 15% vol but the world realizes 30%:
    # a short gamma position bleeds, so mean P&L is negative.
    res = simulate_delta_hedge(100, 100, 1.0, 0.03, 0.15, OptionType.CALL,
                               n_steps=50, n_paths=20_000,
                               hedge_vol=0.15, real_vol=0.30, seed=11)
    assert res.mean_pnl < 0


def test_hedging_at_high_vol_profits_when_realized_is_low():
    # Sell rich (30% priced) and realize cheap (15%): short gamma wins.
    res = simulate_delta_hedge(100, 100, 1.0, 0.03, 0.30, OptionType.CALL,
                               n_steps=50, n_paths=20_000,
                               hedge_vol=0.30, real_vol=0.15, seed=13)
    assert res.mean_pnl > 0


def test_return_samples_length():
    res, samples = simulate_delta_hedge(100, 100, 0.5, 0.05, 0.2,
                                        n_steps=20, n_paths=1000, seed=2,
                                        return_samples=True)
    assert len(samples) == 1000
    assert res.min_pnl <= res.mean_pnl <= res.max_pnl


def test_reproducible_with_seed():
    a = simulate_delta_hedge(100, 100, 1.0, 0.05, 0.2, n_steps=20,
                             n_paths=2000, seed=99)
    b = simulate_delta_hedge(100, 100, 1.0, 0.05, 0.2, n_steps=20,
                             n_paths=2000, seed=99)
    assert a.mean_pnl == b.mean_pnl
    assert a.std_pnl == b.std_pnl
