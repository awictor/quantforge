"""RQMC arithmetic Asian with geometric control variate (sobol_arithmetic_asian_rqmc)."""

import pytest

from quantforge import (
    sobol_arithmetic_asian_rqmc,
    arithmetic_asian_mc,
    OptionType,
)


S, K, T, R, SIG, NS = 100.0, 100.0, 1.0, 0.05, 0.2, 6


@pytest.mark.slow
def test_matches_control_variate_mc():
    cv = sobol_arithmetic_asian_rqmc(S, K, T, R, SIG, OptionType.CALL,
                                     n_steps=NS, n_paths=4096, n_rand=24, seed=1)
    ref = arithmetic_asian_mc(S, K, T, R, SIG, OptionType.CALL, n_steps=NS,
                              n_paths=400_000, control_variate=True, seed=9)
    assert cv.price == pytest.approx(ref.price,
                                     abs=3.0 * (cv.std_error + ref.std_error))


def test_control_variate_lowers_se():
    cv = sobol_arithmetic_asian_rqmc(S, K, T, R, SIG, OptionType.CALL,
                                     control_variate=True, n_steps=NS,
                                     n_paths=2048, n_rand=24, seed=2)
    plain = sobol_arithmetic_asian_rqmc(S, K, T, R, SIG, OptionType.CALL,
                                        control_variate=False, n_steps=NS,
                                        n_paths=2048, n_rand=24, seed=2)
    # The control variate cuts the SE well below the RQMC-only estimator.
    assert cv.std_error < 0.4 * plain.std_error


def test_put_matches_control_variate_mc():
    cv = sobol_arithmetic_asian_rqmc(S, K, T, R, SIG, OptionType.PUT,
                                     n_steps=NS, n_paths=4096, n_rand=24, seed=3)
    ref = arithmetic_asian_mc(S, K, T, R, SIG, OptionType.PUT, n_steps=NS,
                              n_paths=200_000, control_variate=True, seed=9)
    assert cv.price == pytest.approx(ref.price,
                                     abs=3.0 * (cv.std_error + ref.std_error))


def test_n_paths_is_total_points():
    cv = sobol_arithmetic_asian_rqmc(S, K, T, R, SIG, n_steps=4, n_paths=2048,
                                     n_rand=16, seed=4)
    assert cv.n_paths == 2048 * 16


def test_bad_n_rand_raises():
    with pytest.raises(ValueError):
        sobol_arithmetic_asian_rqmc(S, K, T, R, SIG, n_rand=1)


def test_bad_n_steps_raises():
    with pytest.raises(ValueError):
        sobol_arithmetic_asian_rqmc(S, K, T, R, SIG, n_steps=99)
