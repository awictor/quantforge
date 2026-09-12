"""Shout option on a binomial tree."""

import pytest

from quantforge import shout_call, ladder_call
from quantforge.bsm import call_price


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.25


def test_shout_at_least_vanilla():
    assert shout_call(S, K, T, R, SIG, 300) >= call_price(S, K, T, R, SIG)


def test_itm_shout_premium_larger():
    itm_shout = shout_call(130, 100, T, R, SIG, 300)
    itm_van = call_price(130, 100, T, R, SIG)
    assert itm_shout - itm_van > 0


def test_higher_vol_raises_value():
    assert shout_call(S, K, T, R, 0.40, 300) > shout_call(S, K, T, R, SIG, 300)


def test_convergence_across_steps():
    assert abs(shout_call(S, K, T, R, SIG, 100) - shout_call(S, K, T, R, SIG, 300)) < 0.5


def test_validation():
    with pytest.raises(ValueError):
        shout_call(-1, K, T, R, SIG)
    with pytest.raises(ValueError):
        shout_call(S, K, 0, R, SIG)


def test_ladder_no_rungs_equals_tree_vanilla():
    # No rungs -> plain CRR European call. Compare to the same tree, not BSM
    # (the BSM gap is pure discretization error).
    from quantforge import american_price

    tree_van = american_price(S, K, T, R, SIG, option_type="call", b=R, steps=200)
    assert abs(ladder_call(S, K, [], T, R, SIG, 200) - tree_van) < 1e-9


def test_ladder_rungs_below_strike_ignored():
    van = ladder_call(S, K, [], T, R, SIG, 200)
    assert abs(ladder_call(S, K, [80, 90], T, R, SIG, 200) - van) < 1e-9


def test_ladder_at_least_vanilla():
    van = ladder_call(S, K, [], T, R, SIG, 200)
    assert ladder_call(S, K, [110], T, R, SIG, 200) >= van


def test_more_rungs_raise_value():
    one = ladder_call(S, K, [110], T, R, SIG, 200)
    three = ladder_call(S, K, [110, 120, 130], T, R, SIG, 200)
    assert three >= one


def test_ladder_higher_vol_raises_value():
    lo = ladder_call(S, K, [110, 120], T, R, SIG, 200)
    hi = ladder_call(S, K, [110, 120], T, R, 0.40, 200)
    assert hi > lo


def test_ladder_validation():
    with pytest.raises(ValueError):
        ladder_call(-1, K, [110], T, R, SIG)
