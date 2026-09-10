"""RQMC Parisian barrier option (sobol_parisian_rqmc)."""

import pytest

from quantforge import (
    sobol_parisian_rqmc,
    parisian_barrier_mc,
    call_price,
    OptionType,
)


S, K, H, T, R, SIG = 100.0, 100.0, 90.0, 1.0, 0.05, 0.2
WIN, NS = 0.25, 12


@pytest.mark.slow
@pytest.mark.parametrize("barrier", ["down-out", "down-in"])
def test_matches_parisian_barrier_mc(barrier):
    rq = sobol_parisian_rqmc(S, K, H, T, R, SIG, WIN, OptionType.CALL, barrier,
                             n_steps=NS, n_paths=4096, n_rand=24, seed=1)
    ref = parisian_barrier_mc(S, K, H, T, R, SIG, WIN, OptionType.CALL, barrier,
                              n_steps=NS, n_paths=400_000, seed=9)
    assert rq.price == pytest.approx(ref.price,
                                     abs=3.0 * (rq.std_error + ref.std_error))


def test_in_plus_out_equals_vanilla():
    # Parisian knock-in + knock-out partition every path (same activation flag),
    # so at the same seed they sum to the vanilla.
    di = sobol_parisian_rqmc(S, K, H, T, R, SIG, WIN, OptionType.CALL,
                             "down-in", n_steps=NS, n_paths=4096, n_rand=24,
                             seed=2)
    do = sobol_parisian_rqmc(S, K, H, T, R, SIG, WIN, OptionType.CALL,
                             "down-out", n_steps=NS, n_paths=4096, n_rand=24,
                             seed=2)
    assert di.price + do.price == pytest.approx(call_price(S, K, T, R, SIG),
                                                abs=0.05)


def test_longer_window_raises_knockout_value():
    # A longer required window makes the knock-out harder to trigger -> the
    # option survives more often -> higher value.
    short = sobol_parisian_rqmc(S, K, H, T, R, SIG, 0.08, OptionType.CALL,
                                "down-out", n_steps=NS, n_paths=4096, n_rand=16,
                                seed=3).price
    long = sobol_parisian_rqmc(S, K, H, T, R, SIG, 0.5, OptionType.CALL,
                               "down-out", n_steps=NS, n_paths=4096, n_rand=16,
                               seed=3).price
    assert long > short


def test_bad_window_raises():
    with pytest.raises(ValueError):
        sobol_parisian_rqmc(S, K, H, T, R, SIG, 0.0)
    with pytest.raises(ValueError):
        sobol_parisian_rqmc(S, K, H, T, R, SIG, 2.0)


def test_bad_barrier_raises():
    with pytest.raises(ValueError):
        sobol_parisian_rqmc(S, K, H, T, R, SIG, WIN, barrier="sideways")


def test_n_steps_capped_at_twelve():
    with pytest.raises(ValueError):
        sobol_parisian_rqmc(S, K, H, T, R, SIG, WIN, n_steps=13)
