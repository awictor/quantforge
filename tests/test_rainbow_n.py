"""n-asset rainbow options (best-of / worst-of) by Monte Carlo."""

import pytest

from quantforge import rainbow_option_mc as rb
from quantforge.multiasset import best_of_call_closed, worst_of_call_closed
from quantforge.bsm import call_price


def test_two_asset_matches_closed_form():
    args = (100, 100, 100, 1.0, 0.05, 0.2, 0.25, 0.4)
    best = rb([100, 100], 100, 1.0, 0.05, [0.2, 0.25], [[1, 0.4], [0.4, 1]],
              best=True, n_paths=300000, seed=1)
    worst = rb([100, 100], 100, 1.0, 0.05, [0.2, 0.25], [[1, 0.4], [0.4, 1]],
               best=False, n_paths=300000, seed=1)
    assert abs(best - best_of_call_closed(*args)) < 0.2
    assert abs(worst - worst_of_call_closed(*args)) < 0.2


def test_worst_le_single_le_best():
    corr = [[1, 0.4], [0.4, 1]]
    best = rb([100, 100], 100, 1.0, 0.05, [0.2, 0.25], corr, best=True,
              n_paths=300000, seed=1)
    worst = rb([100, 100], 100, 1.0, 0.05, [0.2, 0.25], corr, best=False,
               n_paths=300000, seed=1)
    single = call_price(100, 100, 1.0, 0.05, 0.2)
    assert worst < single < best


def test_three_asset_worst_below_best():
    spots, sig = [100, 95, 105], [0.2, 0.25, 0.3]
    corr = [[1, 0.5, 0.3], [0.5, 1, 0.4], [0.3, 0.4, 1]]
    b = rb(spots, 100, 1.0, 0.05, sig, corr, best=True, n_paths=200000, seed=2)
    w = rb(spots, 100, 1.0, 0.05, sig, corr, best=False, n_paths=200000, seed=2)
    assert w < b


def test_validation():
    with pytest.raises(ValueError):
        rb([100], 100, 1.0, 0.05, [0.2, 0.2], [[1]])       # length mismatch
    with pytest.raises(ValueError):
        rb([100, 100], 100, 1.0, 0.05, [0.2, 0.2], [[1.0]])  # bad corr shape
