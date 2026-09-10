"""Randomized-QMC European price with honest SE (sobol_european_rqmc)."""

import pytest

from quantforge import (
    OptionType,
    sobol_european_rqmc,
    european_mc,
    call_price,
    put_price,
)


S, K, T, R, SIG = 100.0, 105.0, 1.0, 0.04, 0.25


def test_call_matches_black_scholes():
    bs = call_price(S, K, T, R, SIG)
    rq = sobol_european_rqmc(S, K, T, R, SIG, OptionType.CALL,
                             n_paths=4096, n_rand=24, seed=1)
    assert rq.price == pytest.approx(bs, abs=3.0 * rq.std_error)


def test_put_matches_black_scholes():
    bs = put_price(S, K, T, R, SIG)
    rq = sobol_european_rqmc(S, K, T, R, SIG, OptionType.PUT,
                             n_paths=4096, n_rand=24, seed=2)
    assert rq.price == pytest.approx(bs, abs=3.0 * rq.std_error)


def test_n_paths_is_total_points():
    rq = sobol_european_rqmc(S, K, T, R, SIG, n_paths=2048, n_rand=16, seed=3)
    assert rq.n_paths == 2048 * 16


@pytest.mark.slow
def test_rqmc_se_beats_plain_at_equal_points():
    rq = sobol_european_rqmc(S, K, T, R, SIG, OptionType.CALL,
                             n_paths=4096, n_rand=24, seed=4)
    plain = european_mc(S, K, T, R, SIG, OptionType.CALL,
                        n_paths=24 * 4096, seed=4)
    assert rq.std_error < 0.25 * plain.std_error


def test_dividend_carry():
    b = R - 0.03
    bs = call_price(S, K, T, R, SIG, b=b)
    rq = sobol_european_rqmc(S, K, T, R, SIG, OptionType.CALL, b=b,
                             n_paths=4096, n_rand=24, seed=5)
    assert rq.price == pytest.approx(bs, abs=3.0 * rq.std_error)


def test_bad_n_rand_raises():
    with pytest.raises(ValueError):
        sobol_european_rqmc(S, K, T, R, SIG, n_rand=1)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        sobol_european_rqmc(-1, K, T, R, SIG)
