"""RQMC autocallable structured note (sobol_autocallable_rqmc)."""

import pytest

from quantforge import sobol_autocallable_rqmc, autocallable_mc


S, T, R, SIG = 100.0, 3.0, 0.03, 0.25
OBS = [1.0, 2.0, 3.0]
AB, CP, PB = 105.0, 0.07, 70.0


@pytest.mark.slow
def test_matches_autocallable_mc():
    rq = sobol_autocallable_rqmc(S, T, R, SIG, OBS, AB, CP, PB,
                                 n_paths=4096, n_rand=24, seed=1)
    ref = autocallable_mc(S, T, R, SIG, OBS, AB, CP, PB,
                          n_paths=400_000, seed=9)
    assert rq.price == pytest.approx(ref.price,
                                     abs=3.0 * (rq.std_error + ref.std_error))


@pytest.mark.slow
def test_matches_without_protection():
    rq = sobol_autocallable_rqmc(S, T, R, SIG, OBS, AB, CP, None,
                                 n_paths=4096, n_rand=24, seed=2)
    ref = autocallable_mc(S, T, R, SIG, OBS, AB, CP, None,
                          n_paths=400_000, seed=9)
    assert rq.price == pytest.approx(ref.price,
                                     abs=3.0 * (rq.std_error + ref.std_error))


def test_price_positive_and_honest_se_small():
    rq = sobol_autocallable_rqmc(S, T, R, SIG, OBS, AB, CP, PB,
                                 n_paths=4096, n_rand=16, seed=3)
    assert rq.price > 0.0
    assert rq.std_error < 0.01


def test_higher_coupon_raises_value():
    lo = sobol_autocallable_rqmc(S, T, R, SIG, OBS, AB, 0.05, PB,
                                 n_paths=4096, n_rand=16, seed=4).price
    hi = sobol_autocallable_rqmc(S, T, R, SIG, OBS, AB, 0.10, PB,
                                 n_paths=4096, n_rand=16, seed=4).price
    assert hi > lo


def test_n_paths_is_total_points():
    rq = sobol_autocallable_rqmc(S, T, R, SIG, OBS, AB, CP, PB,
                                 n_paths=2048, n_rand=16, seed=5)
    assert rq.n_paths == 2048 * 16


def test_bad_observation_times_raise():
    with pytest.raises(ValueError):
        sobol_autocallable_rqmc(S, T, R, SIG, [1.0, 1.0, 3.0], AB, CP)


def test_last_obs_must_be_maturity():
    with pytest.raises(ValueError):
        sobol_autocallable_rqmc(S, T, R, SIG, [1.0, 2.0], AB, CP)


def test_bad_n_rand_raises():
    with pytest.raises(ValueError):
        sobol_autocallable_rqmc(S, T, R, SIG, OBS, AB, CP, n_rand=1)
