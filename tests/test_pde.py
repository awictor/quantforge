"""Crank-Nicolson PDE solver for European and American options."""

import pytest

from quantforge import (
    OptionType,
    crank_nicolson_price,
    call_price,
    put_price,
    american_price as crr,
)


S, K, T, R, SIGMA = 100.0, 100.0, 1.0, 0.05, 0.2


def test_european_call_matches_black_scholes():
    cn = crank_nicolson_price(S, K, T, R, SIGMA, OptionType.CALL,
                              n_space=300, n_time=300)
    assert cn == pytest.approx(call_price(S, K, T, R, SIGMA), abs=1e-2)


def test_european_put_matches_black_scholes():
    cn = crank_nicolson_price(S, K, T, R, SIGMA, OptionType.PUT,
                              n_space=300, n_time=300)
    assert cn == pytest.approx(put_price(S, K, T, R, SIGMA), abs=1e-2)


def test_convergence_in_grid():
    bs = call_price(S, K, T, R, SIGMA)
    e_lo = abs(crank_nicolson_price(S, K, T, R, SIGMA, n_space=100, n_time=100) - bs)
    e_hi = abs(crank_nicolson_price(S, K, T, R, SIGMA, n_space=400, n_time=400) - bs)
    assert e_hi < e_lo


def test_dividend_carry_matches_black_scholes():
    q = 0.03
    cn = crank_nicolson_price(S, K, T, R, SIGMA, OptionType.CALL, b=R - q,
                              n_space=300, n_time=300)
    assert cn == pytest.approx(call_price(S, K, T, R, SIGMA, b=R - q), abs=1e-2)


def test_flat_local_vol_matches_constant():
    cn = crank_nicolson_price(S, K, T, R, option_type=OptionType.CALL,
                              local_vol_fn=lambda s, t: SIGMA,
                              n_space=300, n_time=300)
    assert cn == pytest.approx(call_price(S, K, T, R, SIGMA), abs=1e-2)


@pytest.mark.slow
def test_american_put_matches_binomial():
    cn = crank_nicolson_price(S, K, T, R, SIGMA, OptionType.PUT, american=True,
                              n_space=300, n_time=300)
    tree = crr(S, K, T, R, SIGMA, OptionType.PUT, b=R, steps=3000)
    assert cn == pytest.approx(tree, abs=1e-2)


def test_american_put_above_european():
    am = crank_nicolson_price(S, K, T, R, SIGMA, OptionType.PUT, american=True,
                              n_space=200, n_time=200)
    eu = crank_nicolson_price(S, K, T, R, SIGMA, OptionType.PUT,
                              n_space=200, n_time=200)
    assert am > eu


def test_zero_time_is_intrinsic():
    assert crank_nicolson_price(120, 100, 0.0, R, SIGMA,
                                OptionType.CALL) == pytest.approx(20.0)


def test_requires_a_vol():
    with pytest.raises(ValueError):
        crank_nicolson_price(S, K, T, R)
