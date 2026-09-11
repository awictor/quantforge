"""RQMC barrier-contingent digital (sobol_barrier_digital_rqmc)."""

import pytest

from quantforge import (
    sobol_barrier_digital_rqmc,
    barrier_digital_mc,
    cash_or_nothing,
    OptionType,
)


S, K, T, R, SIG, NS = 100.0, 100.0, 1.0, 0.05, 0.2, 6


@pytest.mark.slow
@pytest.mark.parametrize("barrier,H", [("up-in", 115.0), ("up-out", 115.0),
                                       ("down-in", 90.0), ("down-out", 90.0)])
def test_matches_barrier_digital_mc(barrier, H):
    rq = sobol_barrier_digital_rqmc(S, K, H, T, R, SIG, OptionType.CALL, barrier,
                                    n_steps=NS, n_paths=4096, n_rand=24, seed=1)
    ref = barrier_digital_mc(S, K, H, T, R, SIG, OptionType.CALL, barrier,
                             n_steps=NS, n_paths=300_000, seed=9)
    assert rq.price == pytest.approx(ref.price,
                                     abs=3.0 * (rq.std_error + ref.std_error))


def test_in_plus_out_equals_unconditional_digital():
    # Same seed -> same paths; knock-in + knock-out partition every path, so
    # they sum to the plain cash-or-nothing digital.
    # Partition identity is exact at a shared seed for any n_rand -- keep small.
    di = sobol_barrier_digital_rqmc(S, K, 90.0, T, R, SIG, OptionType.CALL,
                                    "down-in", n_steps=NS, n_paths=4096,
                                    n_rand=6, seed=2)
    do = sobol_barrier_digital_rqmc(S, K, 90.0, T, R, SIG, OptionType.CALL,
                                    "down-out", n_steps=NS, n_paths=4096,
                                    n_rand=6, seed=2)
    dig = cash_or_nothing(S, K, T, R, SIG, OptionType.CALL)
    assert di.price + do.price == pytest.approx(dig, abs=5e-3)


def test_cash_scales_linearly():
    a = sobol_barrier_digital_rqmc(S, K, 90.0, T, R, SIG, OptionType.CALL,
                                   "down-out", cash=1.0, n_steps=NS,
                                   n_paths=4096, n_rand=12, seed=3)
    b = sobol_barrier_digital_rqmc(S, K, 90.0, T, R, SIG, OptionType.CALL,
                                   "down-out", cash=8.0, n_steps=NS,
                                   n_paths=4096, n_rand=12, seed=3)
    assert b.price == pytest.approx(8.0 * a.price, rel=1e-9)


def test_n_paths_is_total_points():
    rq = sobol_barrier_digital_rqmc(S, K, 90.0, T, R, SIG, barrier="down-out",
                                    n_steps=4, n_paths=2048, n_rand=16, seed=4)
    assert rq.n_paths == 2048 * 16


def test_bad_barrier_raises():
    with pytest.raises(ValueError):
        sobol_barrier_digital_rqmc(S, K, 90.0, T, R, SIG, barrier="sideways")


def test_bad_n_rand_raises():
    with pytest.raises(ValueError):
        sobol_barrier_digital_rqmc(S, K, 90.0, T, R, SIG, n_rand=1)
