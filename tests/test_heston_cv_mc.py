"""Heston QE Monte Carlo with underlying control variate (heston_cv_mc)."""

import pytest

from quantforge import (
    OptionType,
    heston_cv_mc,
    heston_qe_mc,
    heston_price,
)


S, K, T, R = 100.0, 100.0, 1.0, 0.03
V0, KAPPA, THETA, XI, RHO = 0.04, 1.5, 0.04, 0.5, -0.7


@pytest.mark.slow
def test_call_matches_fourier():
    fx = heston_price(S, K, T, R, V0, KAPPA, THETA, XI, RHO, OptionType.CALL)
    cv = heston_cv_mc(S, K, T, R, V0, KAPPA, THETA, XI, RHO, OptionType.CALL,
                      n_paths=120_000, seed=1)
    assert cv.price == pytest.approx(fx, abs=3.0 * cv.std_error)


@pytest.mark.slow
def test_put_matches_fourier():
    fx = heston_price(S, K, T, R, V0, KAPPA, THETA, XI, RHO, OptionType.PUT)
    cv = heston_cv_mc(S, K, T, R, V0, KAPPA, THETA, XI, RHO, OptionType.PUT,
                      n_paths=120_000, seed=2)
    assert cv.price == pytest.approx(fx, abs=3.0 * cv.std_error)


def test_control_variate_lowers_standard_error():
    cv = heston_cv_mc(S, K, T, R, V0, KAPPA, THETA, XI, RHO, OptionType.CALL,
                      n_paths=30_000, seed=3)
    plain = heston_qe_mc(S, K, T, R, V0, KAPPA, THETA, XI, RHO, OptionType.CALL,
                         n_paths=30_000, seed=3)
    assert cv.std_error < plain.std_error


def test_dividend_yield():
    q = 0.02
    fx = heston_price(S, K, T, R, V0, KAPPA, THETA, XI, RHO, OptionType.CALL, q=q)
    cv = heston_cv_mc(S, K, T, R, V0, KAPPA, THETA, XI, RHO, OptionType.CALL,
                      q=q, n_paths=60_000, seed=4)
    assert cv.price == pytest.approx(fx, abs=3.0 * cv.std_error)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        heston_cv_mc(-1, K, T, R, V0, KAPPA, THETA, XI, RHO)
