"""Randomized-QMC discrete fixed-strike lookback (sobol_fixed_lookback_rqmc)."""

import pytest

from quantforge import (
    sobol_fixed_lookback_rqmc,
    fixed_strike_lookback,
    OptionType,
)


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2


def test_discrete_below_continuous_call():
    cont = fixed_strike_lookback(S, K, T, R, SIG, OptionType.CALL)
    mc = sobol_fixed_lookback_rqmc(S, K, T, R, SIG, OptionType.CALL, n_steps=6,
                                   n_paths=4096, n_rand=24, seed=1)
    assert mc.price < cont
    assert mc.price > 0.0


def test_discrete_below_continuous_put():
    cont = fixed_strike_lookback(S, K, T, R, SIG, OptionType.PUT)
    mc = sobol_fixed_lookback_rqmc(S, K, T, R, SIG, OptionType.PUT, n_steps=6,
                                   n_paths=4096, n_rand=24, seed=2)
    assert mc.price < cont
    assert mc.price > 0.0


def test_price_rises_with_monitoring_frequency():
    p2 = sobol_fixed_lookback_rqmc(S, K, T, R, SIG, OptionType.CALL, n_steps=2,
                                   n_paths=4096, n_rand=16, seed=3).price
    p6 = sobol_fixed_lookback_rqmc(S, K, T, R, SIG, OptionType.CALL, n_steps=6,
                                   n_paths=4096, n_rand=16, seed=3).price
    assert p6 > p2


def test_call_at_least_vanilla():
    # A fixed-strike lookback call is worth at least the plain vanilla (its max
    # over the path dominates the terminal spot).
    from quantforge import call_price
    mc = sobol_fixed_lookback_rqmc(S, K, T, R, SIG, OptionType.CALL, n_steps=6,
                                   n_paths=4096, n_rand=24, seed=4)
    assert mc.price > call_price(S, K, T, R, SIG) - 3 * mc.std_error


def test_n_paths_is_total_points():
    mc = sobol_fixed_lookback_rqmc(S, K, T, R, SIG, n_steps=4, n_paths=2048,
                                   n_rand=16, seed=5)
    assert mc.n_paths == 2048 * 16


def test_bad_n_rand_raises():
    with pytest.raises(ValueError):
        sobol_fixed_lookback_rqmc(S, K, T, R, SIG, n_rand=1)


def test_bad_n_steps_raises():
    with pytest.raises(ValueError):
        sobol_fixed_lookback_rqmc(S, K, T, R, SIG, n_steps=99)
