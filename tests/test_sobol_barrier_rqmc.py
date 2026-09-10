"""Randomized-QMC discretely-monitored single-barrier option (sobol_barrier_rqmc)."""

import pytest

from quantforge import (
    sobol_barrier_rqmc,
    barrier_mc,
    call_price,
    OptionType,
)


S, K, T, R, SIG, NS = 100.0, 100.0, 1.0, 0.05, 0.2, 6


@pytest.mark.slow
@pytest.mark.parametrize("barrier,H", [("down-out", 90.0), ("down-in", 90.0),
                                       ("up-out", 115.0), ("up-in", 115.0)])
def test_matches_discrete_barrier_mc(barrier, H):
    rq = sobol_barrier_rqmc(S, K, H, T, R, SIG, OptionType.CALL, barrier,
                            n_steps=NS, n_paths=4096, n_rand=24, seed=1)
    ref = barrier_mc(S, K, H, T, R, SIG, OptionType.CALL, barrier, n_steps=NS,
                     n_paths=300_000, seed=9, brownian_bridge=False)
    assert rq.price == pytest.approx(ref.price,
                                     abs=3.0 * (rq.std_error + ref.std_error))


def test_in_plus_out_equals_vanilla():
    # Same seed -> same paths; knock-in + knock-out partition every path, so
    # they sum to the (discretely-sampled) vanilla, which equals the analytic
    # European exactly since the terminal spot is the same.
    di = sobol_barrier_rqmc(S, K, 90.0, T, R, SIG, OptionType.CALL, "down-in",
                            n_steps=NS, n_paths=4096, n_rand=24, seed=2)
    do = sobol_barrier_rqmc(S, K, 90.0, T, R, SIG, OptionType.CALL, "down-out",
                            n_steps=NS, n_paths=4096, n_rand=24, seed=2)
    vanilla = call_price(S, K, T, R, SIG)
    assert di.price + do.price == pytest.approx(vanilla, abs=0.05)


def test_knockout_below_vanilla():
    do = sobol_barrier_rqmc(S, K, 90.0, T, R, SIG, OptionType.CALL, "down-out",
                            n_steps=NS, n_paths=4096, n_rand=24, seed=3)
    assert do.price < call_price(S, K, T, R, SIG)


def test_n_paths_is_total_points():
    rq = sobol_barrier_rqmc(S, K, 90.0, T, R, SIG, barrier="down-out",
                            n_steps=4, n_paths=2048, n_rand=16, seed=4)
    assert rq.n_paths == 2048 * 16


def test_bad_barrier_raises():
    with pytest.raises(ValueError):
        sobol_barrier_rqmc(S, K, 90.0, T, R, SIG, barrier="sideways")


def test_bad_H_raises():
    with pytest.raises(ValueError):
        sobol_barrier_rqmc(S, K, -1.0, T, R, SIG)


def test_bad_n_rand_raises():
    with pytest.raises(ValueError):
        sobol_barrier_rqmc(S, K, 90.0, T, R, SIG, n_rand=1)
