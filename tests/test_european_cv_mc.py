"""European MC with underlying control variate + antithetic (european_cv_mc)."""

import pytest

from quantforge import (
    OptionType,
    european_cv_mc,
    european_mc,
    call_price,
    put_price,
)


S, K, T, R, SIG = 100.0, 105.0, 1.0, 0.04, 0.25


def test_call_matches_black_scholes():
    bs = call_price(S, K, T, R, SIG)
    cv = european_cv_mc(S, K, T, R, SIG, OptionType.CALL, n_paths=60_000, seed=1)
    assert cv.price == pytest.approx(bs, abs=3.0 * cv.std_error)


def test_put_matches_black_scholes():
    bs = put_price(S, K, T, R, SIG)
    cv = european_cv_mc(S, K, T, R, SIG, OptionType.PUT, n_paths=60_000, seed=2)
    assert cv.price == pytest.approx(bs, abs=3.0 * cv.std_error)


def test_control_variate_lowers_standard_error():
    # Same path count, same seed: the control variate must reduce the SE.
    cv = european_cv_mc(S, K, T, R, SIG, OptionType.CALL, n_paths=40_000, seed=3)
    plain = european_mc(S, K, T, R, SIG, OptionType.CALL, n_paths=40_000, seed=3)
    assert cv.std_error < plain.std_error


@pytest.mark.slow
def test_control_variate_beats_antithetic_only():
    # At-the-money, where correlation between payoff and terminal spot is high,
    # the SE reduction should be substantial (well under half the plain SE).
    cv = european_cv_mc(100.0, 100.0, T, R, SIG, n_paths=80_000, seed=4)
    plain = european_mc(100.0, 100.0, T, R, SIG, n_paths=80_000, seed=4)
    assert cv.std_error < 0.5 * plain.std_error


def test_dividend_yield_carry():
    b = R - 0.03                       # cost of carry with a 3% dividend yield
    bs = call_price(S, K, T, R, SIG, b=b)
    cv = european_cv_mc(S, K, T, R, SIG, OptionType.CALL, b=b,
                        n_paths=60_000, seed=5)
    assert cv.price == pytest.approx(bs, abs=3.0 * cv.std_error)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        european_cv_mc(-1, K, T, R, SIG)
