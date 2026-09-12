"""Shout option on a binomial tree."""

import pytest

from quantforge import shout_call
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
