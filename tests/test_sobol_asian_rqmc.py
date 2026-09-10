"""Randomized-QMC arithmetic Asian with honest SE (sobol_asian_rqmc)."""

import pytest

from quantforge import (
    OptionType,
    sobol_asian_rqmc,
    arithmetic_asian_mc,
)


S, K, T, R, SIG, NS = 100.0, 100.0, 1.0, 0.05, 0.2, 6


def test_call_matches_control_variate_mc():
    ref = arithmetic_asian_mc(S, K, T, R, SIG, OptionType.CALL, n_steps=NS,
                              n_paths=120_000, control_variate=True, seed=9)
    rq = sobol_asian_rqmc(S, K, T, R, SIG, OptionType.CALL, n_steps=NS,
                          n_paths=4096, n_rand=16, seed=1)
    assert rq.price == pytest.approx(ref.price,
                                     abs=3.0 * (ref.std_error + rq.std_error))


def test_put_matches_control_variate_mc():
    ref = arithmetic_asian_mc(S, K, T, R, SIG, OptionType.PUT, n_steps=NS,
                              n_paths=120_000, control_variate=True, seed=9)
    rq = sobol_asian_rqmc(S, K, T, R, SIG, OptionType.PUT, n_steps=NS,
                          n_paths=4096, n_rand=16, seed=2)
    assert rq.price == pytest.approx(ref.price,
                                     abs=3.0 * (ref.std_error + rq.std_error))


def test_n_paths_is_total_points():
    rq = sobol_asian_rqmc(S, K, T, R, SIG, n_steps=4, n_paths=2048, n_rand=16,
                          seed=3)
    assert rq.n_paths == 2048 * 16


@pytest.mark.slow
def test_rqmc_asian_se_beats_plain():
    rq = sobol_asian_rqmc(S, K, T, R, SIG, OptionType.CALL, n_steps=NS,
                          n_paths=4096, n_rand=24, seed=4)
    plain = arithmetic_asian_mc(S, K, T, R, SIG, OptionType.CALL, n_steps=NS,
                                n_paths=24 * 4096, control_variate=False, seed=4)
    assert rq.std_error < 0.25 * plain.std_error


def test_bad_n_rand_raises():
    with pytest.raises(ValueError):
        sobol_asian_rqmc(S, K, T, R, SIG, n_rand=1)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        sobol_asian_rqmc(-1, K, T, R, SIG)
