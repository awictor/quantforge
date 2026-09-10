"""European MC by importance sampling for deep-OTM options (european_is_mc)."""

import pytest

from quantforge import (
    OptionType,
    european_is_mc,
    european_mc,
    call_price,
    put_price,
)


S, T, R, SIG = 100.0, 1.0, 0.03, 0.2


def test_deep_otm_call_matches_black_scholes():
    K = 160.0
    bs = call_price(S, K, T, R, SIG)
    mc = european_is_mc(S, K, T, R, SIG, OptionType.CALL, n_paths=60_000, seed=1)
    assert mc.price == pytest.approx(bs, abs=3.0 * mc.std_error)


def test_deep_otm_put_matches_black_scholes():
    K = 55.0
    bs = put_price(S, K, T, R, SIG)
    mc = european_is_mc(S, K, T, R, SIG, OptionType.PUT, n_paths=60_000, seed=2)
    assert mc.price == pytest.approx(bs, abs=3.0 * mc.std_error)


def test_importance_sampling_beats_plain_deep_otm():
    # Deep OTM is where importance sampling shines: a large SE reduction.
    K = 160.0
    is_ = european_is_mc(S, K, T, R, SIG, OptionType.CALL, n_paths=40_000, seed=3)
    plain = european_mc(S, K, T, R, SIG, OptionType.CALL, n_paths=40_000, seed=3)
    assert is_.std_error < 0.25 * plain.std_error


def test_atm_still_correct():
    K = 100.0
    bs = call_price(S, K, T, R, SIG)
    mc = european_is_mc(S, K, T, R, SIG, OptionType.CALL, n_paths=60_000, seed=4)
    assert mc.price == pytest.approx(bs, abs=3.0 * mc.std_error)


def test_explicit_shift_still_unbiased():
    K = 130.0
    bs = call_price(S, K, T, R, SIG)
    mc = european_is_mc(S, K, T, R, SIG, OptionType.CALL, shift=1.0,
                        n_paths=80_000, seed=5)
    assert mc.price == pytest.approx(bs, abs=3.0 * mc.std_error)


def test_dividend_carry():
    K = 150.0
    b = R - 0.04
    bs = call_price(S, K, T, R, SIG, b=b)
    mc = european_is_mc(S, K, T, R, SIG, OptionType.CALL, b=b,
                        n_paths=60_000, seed=6)
    assert mc.price == pytest.approx(bs, abs=3.0 * mc.std_error)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        european_is_mc(-1, 160.0, T, R, SIG)
