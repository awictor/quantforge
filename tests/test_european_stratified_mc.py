"""European MC by stratified sampling of the terminal normal."""

import pytest

from quantforge import (
    OptionType,
    european_stratified_mc,
    european_mc,
    call_price,
    put_price,
)


S, K, T, R, SIG = 100.0, 105.0, 1.0, 0.04, 0.25


def test_call_matches_black_scholes():
    bs = call_price(S, K, T, R, SIG)
    mc = european_stratified_mc(S, K, T, R, SIG, OptionType.CALL,
                                n_strata=400, n_per=25, seed=1)
    assert mc.price == pytest.approx(bs, abs=3.0 * mc.std_error)


def test_put_matches_black_scholes():
    bs = put_price(S, K, T, R, SIG)
    mc = european_stratified_mc(S, K, T, R, SIG, OptionType.PUT,
                                n_strata=400, n_per=25, seed=2)
    assert mc.price == pytest.approx(bs, abs=3.0 * mc.std_error)


def test_stratification_lowers_standard_error():
    # 10k paths each: stratified must beat plain by a wide margin on a smooth
    # payoff.
    st = european_stratified_mc(S, K, T, R, SIG, OptionType.CALL,
                                n_strata=400, n_per=25, seed=3)
    plain = european_mc(S, K, T, R, SIG, OptionType.CALL,
                        n_paths=10_000, seed=3)
    assert st.n_paths == 10_000
    assert st.std_error < 0.25 * plain.std_error


def test_dividend_carry():
    b = R - 0.03
    bs = call_price(S, K, T, R, SIG, b=b)
    mc = european_stratified_mc(S, K, T, R, SIG, OptionType.CALL, b=b,
                                n_strata=400, n_per=25, seed=4)
    assert mc.price == pytest.approx(bs, abs=3.0 * mc.std_error)


def test_deep_otm_still_unbiased():
    K2 = 160.0
    bs = call_price(S, K2, T, R, SIG)
    mc = european_stratified_mc(S, K2, T, R, SIG, OptionType.CALL,
                                n_strata=1000, n_per=20, seed=5)
    assert mc.price == pytest.approx(bs, abs=3.0 * mc.std_error)


def test_bad_n_per_raises():
    with pytest.raises(ValueError):
        european_stratified_mc(S, K, T, R, SIG, n_per=1)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        european_stratified_mc(-1, K, T, R, SIG)
