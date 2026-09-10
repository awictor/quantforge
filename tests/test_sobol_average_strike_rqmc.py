"""RQMC average-strike Asian option (sobol_average_strike_rqmc)."""

import pytest

from quantforge import (
    sobol_average_strike_rqmc,
    average_strike_asian_mc,
    OptionType,
)


S, T, R, SIG, NS = 100.0, 1.0, 0.05, 0.2, 6


@pytest.mark.slow
@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
def test_matches_average_strike_mc(ot):
    rq = sobol_average_strike_rqmc(S, T, R, SIG, ot, n_steps=NS,
                                   n_paths=4096, n_rand=24, seed=1)
    ref = average_strike_asian_mc(S, T, R, SIG, ot, n_steps=NS,
                                  n_paths=400_000, seed=9)
    assert rq.price == pytest.approx(ref.price,
                                     abs=3.0 * (rq.std_error + ref.std_error))


def test_prices_positive():
    c = sobol_average_strike_rqmc(S, T, R, SIG, OptionType.CALL, n_steps=NS,
                                  n_paths=4096, n_rand=16, seed=2)
    p = sobol_average_strike_rqmc(S, T, R, SIG, OptionType.PUT, n_steps=NS,
                                  n_paths=4096, n_rand=16, seed=3)
    assert c.price > 0.0
    assert p.price > 0.0


def test_honest_se_small():
    rq = sobol_average_strike_rqmc(S, T, R, SIG, OptionType.CALL, n_steps=NS,
                                   n_paths=4096, n_rand=24, seed=4)
    assert rq.std_error < 0.02


def test_n_paths_is_total_points():
    rq = sobol_average_strike_rqmc(S, T, R, SIG, n_steps=4, n_paths=2048,
                                   n_rand=16, seed=5)
    assert rq.n_paths == 2048 * 16


def test_bad_n_rand_raises():
    with pytest.raises(ValueError):
        sobol_average_strike_rqmc(S, T, R, SIG, n_rand=1)


def test_bad_n_steps_raises():
    with pytest.raises(ValueError):
        sobol_average_strike_rqmc(S, T, R, SIG, n_steps=99)
