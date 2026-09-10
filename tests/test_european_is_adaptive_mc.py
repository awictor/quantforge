"""Adaptive (pilot-tuned) importance-sampling European MC."""

import pytest

from quantforge import (
    OptionType,
    european_is_adaptive_mc,
    european_is_mc,
    call_price,
    put_price,
)


S, T, R, SIG = 100.0, 1.0, 0.03, 0.2


def test_deep_otm_call_matches_black_scholes():
    K = 160.0
    bs = call_price(S, K, T, R, SIG)
    mc = european_is_adaptive_mc(S, K, T, R, SIG, OptionType.CALL,
                                 n_paths=60_000, seed=1)
    assert mc.price == pytest.approx(bs, abs=3.0 * mc.std_error)


def test_deep_otm_put_matches_black_scholes():
    K = 55.0
    bs = put_price(S, K, T, R, SIG)
    mc = european_is_adaptive_mc(S, K, T, R, SIG, OptionType.PUT,
                                 n_paths=60_000, seed=2)
    assert mc.price == pytest.approx(bs, abs=3.0 * mc.std_error)


@pytest.mark.slow
def test_adaptive_shift_at_or_below_fixed_shift():
    # The pilot-tuned shift should never be worse than the strike-centring one
    # for a deep-OTM vanilla (whose optimal shift sits further OTM).
    K = 160.0
    ad = european_is_adaptive_mc(S, K, T, R, SIG, OptionType.CALL,
                                 n_pilot=40_000, n_paths=80_000, seed=3)
    fx = european_is_mc(S, K, T, R, SIG, OptionType.CALL,
                        n_paths=80_000, seed=3)
    assert ad.std_error <= fx.std_error * 1.02


def test_atm_still_correct():
    K = 100.0
    bs = call_price(S, K, T, R, SIG)
    mc = european_is_adaptive_mc(S, K, T, R, SIG, OptionType.CALL,
                                 n_paths=60_000, seed=4)
    assert mc.price == pytest.approx(bs, abs=3.0 * mc.std_error)


def test_dividend_carry():
    K = 150.0
    b = R - 0.04
    bs = call_price(S, K, T, R, SIG, b=b)
    mc = european_is_adaptive_mc(S, K, T, R, SIG, OptionType.CALL, b=b,
                                 n_paths=60_000, seed=5)
    assert mc.price == pytest.approx(bs, abs=3.0 * mc.std_error)


def test_bad_pilot_raises():
    with pytest.raises(ValueError):
        european_is_adaptive_mc(S, 160.0, T, R, SIG, n_pilot=10)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        european_is_adaptive_mc(-1, 160.0, T, R, SIG)
