"""RQMC double-knockout corridor option (sobol_double_knockout_rqmc)."""

import pytest

from quantforge import (
    sobol_double_knockout_rqmc,
    double_knockout_mc,
    call_price,
    OptionType,
)


S, K, T, R, SIG, NS = 100.0, 100.0, 1.0, 0.05, 0.2, 6
LO, HI = 85.0, 120.0


@pytest.mark.slow
def test_matches_double_knockout_mc():
    rq = sobol_double_knockout_rqmc(S, K, T, R, SIG, LO, HI, OptionType.CALL,
                                    n_steps=NS, n_paths=4096, n_rand=24, seed=1)
    ref = double_knockout_mc(S, K, T, R, SIG, LO, HI, OptionType.CALL,
                             n_steps=NS, n_paths=400_000, seed=9)
    assert rq.price == pytest.approx(ref.price,
                                     abs=3.0 * (rq.std_error + ref.std_error))


def test_wide_corridor_near_vanilla():
    # A very wide corridor is rarely breached -> close to the vanilla call.
    w = sobol_double_knockout_rqmc(S, K, T, R, SIG, 40.0, 250.0, OptionType.CALL,
                                   n_steps=NS, n_paths=4096, n_rand=24, seed=2)
    assert w.price == pytest.approx(call_price(S, K, T, R, SIG), abs=0.2)


def test_narrower_corridor_worth_less():
    narrow = sobol_double_knockout_rqmc(S, K, T, R, SIG, LO, HI, OptionType.CALL,
                                        n_steps=NS, n_paths=4096, n_rand=16,
                                        seed=3).price
    wide = sobol_double_knockout_rqmc(S, K, T, R, SIG, 60.0, 160.0,
                                      OptionType.CALL, n_steps=NS, n_paths=4096,
                                      n_rand=16, seed=3).price
    assert narrow < wide


def test_n_paths_is_total_points():
    rq = sobol_double_knockout_rqmc(S, K, T, R, SIG, LO, HI, n_steps=4,
                                    n_paths=2048, n_rand=16, seed=4)
    assert rq.n_paths == 2048 * 16


def test_spot_outside_corridor_raises():
    with pytest.raises(ValueError):
        sobol_double_knockout_rqmc(S, K, T, R, SIG, 105.0, 120.0)


def test_bad_n_rand_raises():
    with pytest.raises(ValueError):
        sobol_double_knockout_rqmc(S, K, T, R, SIG, LO, HI, n_rand=1)
