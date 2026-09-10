"""RQMC capped cliquet (ratchet) note (sobol_cliquet_rqmc)."""

import pytest

from quantforge import sobol_cliquet_rqmc, capped_cliquet_mc


S, T, R, SIG = 100.0, 1.0, 0.05, 0.25
RT = [0.25, 0.5, 0.75, 1.0]
LC, LF, GC, GF = 0.05, -0.03, 0.20, 0.0


@pytest.mark.slow
def test_matches_capped_cliquet_mc():
    rq = sobol_cliquet_rqmc(S, T, R, SIG, RT, LC, LF, GC, GF,
                            n_paths=4096, n_rand=24, seed=1)
    ref = capped_cliquet_mc(S, T, R, SIG, RT, LC, LF, GC, GF,
                            n_paths=400_000, seed=9)
    assert rq.price == pytest.approx(ref.price,
                                     abs=3.0 * (rq.std_error + ref.std_error))


@pytest.mark.slow
def test_matches_uncapped_local_no_global():
    rq = sobol_cliquet_rqmc(S, T, R, SIG, RT, None, LF, None, GF,
                            n_paths=4096, n_rand=24, seed=2)
    ref = capped_cliquet_mc(S, T, R, SIG, RT, None, LF, None, GF,
                            n_paths=400_000, seed=9)
    assert rq.price == pytest.approx(ref.price,
                                     abs=3.0 * (rq.std_error + ref.std_error))


def test_local_cap_lowers_value():
    tight = sobol_cliquet_rqmc(S, T, R, SIG, RT, 0.02, LF, GC, GF,
                               n_paths=4096, n_rand=16, seed=3).price
    loose = sobol_cliquet_rqmc(S, T, R, SIG, RT, 0.10, LF, GC, GF,
                               n_paths=4096, n_rand=16, seed=3).price
    assert tight < loose


def test_price_positive_and_honest_se_small():
    rq = sobol_cliquet_rqmc(S, T, R, SIG, RT, LC, LF, GC, GF,
                            n_paths=4096, n_rand=16, seed=4)
    assert rq.price > 0.0
    assert rq.std_error < 0.01


def test_n_paths_is_total_points():
    rq = sobol_cliquet_rqmc(S, T, R, SIG, RT, LC, LF, GC, GF,
                            n_paths=2048, n_rand=16, seed=5)
    assert rq.n_paths == 2048 * 16


def test_bad_reset_times_raise():
    with pytest.raises(ValueError):
        sobol_cliquet_rqmc(S, T, R, SIG, [0.5, 0.5, 1.0])


def test_last_reset_must_be_maturity():
    with pytest.raises(ValueError):
        sobol_cliquet_rqmc(S, T, R, SIG, [0.25, 0.5])


def test_bad_n_rand_raises():
    with pytest.raises(ValueError):
        sobol_cliquet_rqmc(S, T, R, SIG, RT, n_rand=1)
