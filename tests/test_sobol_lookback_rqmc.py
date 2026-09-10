"""Randomized-QMC discrete floating-strike lookback (sobol_lookback_rqmc)."""

import pytest

from quantforge import (
    sobol_lookback_rqmc,
    floating_strike_lookback,
    OptionType,
)


S, T, R, SIG = 100.0, 1.0, 0.05, 0.2


def test_discrete_below_continuous_call():
    # Discrete monitoring sees less extreme lows -> lookback worth less than the
    # continuously-monitored Goldman-Sosin-Gatto closed form.
    cont = floating_strike_lookback(S, T, R, SIG, OptionType.CALL)
    mc = sobol_lookback_rqmc(S, T, R, SIG, OptionType.CALL, n_steps=6,
                             n_paths=2048, n_rand=12, seed=1)
    assert mc.price < cont
    assert mc.price > 0.0


def test_discrete_below_continuous_put():
    cont = floating_strike_lookback(S, T, R, SIG, OptionType.PUT)
    mc = sobol_lookback_rqmc(S, T, R, SIG, OptionType.PUT, n_steps=6,
                             n_paths=2048, n_rand=12, seed=2)
    assert mc.price < cont
    assert mc.price > 0.0


def test_price_rises_with_monitoring_frequency():
    # More monitoring dates capture more extreme lows -> higher call value,
    # approaching the continuous price.
    p2 = sobol_lookback_rqmc(S, T, R, SIG, OptionType.CALL, n_steps=2,
                             n_paths=4096, n_rand=16, seed=3).price
    p6 = sobol_lookback_rqmc(S, T, R, SIG, OptionType.CALL, n_steps=6,
                             n_paths=4096, n_rand=16, seed=3).price
    assert p6 > p2


def test_n_paths_is_total_points():
    mc = sobol_lookback_rqmc(S, T, R, SIG, n_steps=4, n_paths=2048, n_rand=16,
                             seed=4)
    assert mc.n_paths == 2048 * 16


def test_honest_se_small():
    # Bridge + RQMC gives a tight across-randomization SE for this smooth-ish
    # path functional.
    mc = sobol_lookback_rqmc(S, T, R, SIG, OptionType.CALL, n_steps=6,
                             n_paths=2048, n_rand=12, seed=5)
    assert mc.std_error < 0.05


def test_bad_n_rand_raises():
    with pytest.raises(ValueError):
        sobol_lookback_rqmc(S, T, R, SIG, n_rand=1)


def test_bad_n_steps_raises():
    with pytest.raises(ValueError):
        sobol_lookback_rqmc(S, T, R, SIG, n_steps=99)
