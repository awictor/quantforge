"""Multi-level Monte Carlo for the arithmetic-average Asian option."""

import pytest

from quantforge import OptionType, mlmc_asian, arithmetic_asian


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2


@pytest.mark.slow
def test_converges_to_turnbull_wakeman():
    tw = arithmetic_asian(S, K, T, R, SIG, OptionType.CALL)
    ml = mlmc_asian(S, K, T, R, SIG, OptionType.CALL, levels=6, n_paths=80_000,
                    seed=3)
    assert ml.price == pytest.approx(tw, abs=3.0 * ml.std_error + 0.05)


@pytest.mark.slow
def test_more_levels_reduce_discretisation_bias():
    tw = arithmetic_asian(S, K, T, R, SIG, OptionType.CALL)
    coarse = mlmc_asian(S, K, T, R, SIG, OptionType.CALL, levels=2,
                        n_paths=80_000, seed=5)
    fine = mlmc_asian(S, K, T, R, SIG, OptionType.CALL, levels=6,
                      n_paths=80_000, seed=5)
    assert abs(fine.price - tw) < abs(coarse.price - tw)


def test_positive_price():
    ml = mlmc_asian(S, K, T, R, SIG, OptionType.CALL, levels=4, n_paths=20_000,
                    seed=1)
    assert ml.price > 0.0
    assert ml.n_paths > 20_000     # levels add samples


def test_put_positive():
    ml = mlmc_asian(S, K, T, R, SIG, OptionType.PUT, levels=4, n_paths=20_000,
                    seed=2)
    assert ml.price > 0.0


def test_bad_params_raise():
    with pytest.raises(ValueError):
        mlmc_asian(S, K, T, R, SIG, levels=-1)
    with pytest.raises(ValueError):
        mlmc_asian(S, K, T, R, SIG, M=1)
