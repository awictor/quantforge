"""Crank-Nicolson barrier PDE (absorbing boundary) vs the closed form."""

import pytest

from quantforge import (
    OptionType,
    crank_nicolson_barrier,
    barrier_option,
    call_price,
)
from quantforge.exotics import Barrier


S, K, T, R, SIGMA = 100.0, 100.0, 1.0, 0.05, 0.25


CASES = [
    ("down-out", Barrier.DOWN_OUT, 90, OptionType.CALL),
    ("down-in", Barrier.DOWN_IN, 90, OptionType.CALL),
    ("up-out", Barrier.UP_OUT, 130, OptionType.CALL),
    ("up-in", Barrier.UP_IN, 130, OptionType.CALL),
    ("down-out", Barrier.DOWN_OUT, 90, OptionType.PUT),
]


@pytest.mark.slow
@pytest.mark.parametrize("name,kind,H,ot", CASES)
def test_barrier_pde_close_to_closed_form(name, kind, H, ot):
    pde = crank_nicolson_barrier(S, K, H, T, R, SIGMA, ot, name,
                                 n_space=1000, n_time=1000)
    cf = barrier_option(S, K, H, T, R, SIGMA, ot, kind, b=R)
    # Node-absorbing barrier is O(ds); allow the known discretisation bias.
    assert pde == pytest.approx(cf, abs=0.15)


def test_in_out_parity_is_vanilla():
    ki = crank_nicolson_barrier(S, K, 90, T, R, SIGMA, OptionType.CALL,
                                "down-in", n_space=600, n_time=600)
    ko = crank_nicolson_barrier(S, K, 90, T, R, SIGMA, OptionType.CALL,
                                "down-out", n_space=600, n_time=600)
    van = call_price(S, K, T, R, SIGMA)
    assert (ki + ko) == pytest.approx(van, abs=1e-2)


def test_knockout_below_vanilla():
    ko = crank_nicolson_barrier(S, K, 90, T, R, SIGMA, OptionType.CALL,
                                "down-out", n_space=400, n_time=400)
    van = call_price(S, K, T, R, SIGMA)
    assert 0.0 < ko < van


def test_spot_at_barrier_is_rebate():
    # A down-and-out with spot already at/below the barrier is dead -> rebate.
    v = crank_nicolson_barrier(90, K, 90, T, R, SIGMA, OptionType.CALL,
                               "down-out", rebate=0.0, n_space=200, n_time=200)
    assert v == pytest.approx(0.0, abs=0.15)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        crank_nicolson_barrier(S, K, -1, T, R, SIGMA, barrier="down-out")
    with pytest.raises(ValueError):
        crank_nicolson_barrier(S, K, 90, T, R, barrier="down-out")  # no vol
    with pytest.raises(ValueError):
        crank_nicolson_barrier(S, K, 90, T, R, SIGMA, barrier="sideways")
