"""Mixed pathwise-likelihood-ratio gamma vs Black-Scholes."""

import pytest

from quantforge import mixed_gamma, lr_greeks, gamma as bs_gamma, OptionType


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2


@pytest.mark.slow
def test_matches_black_scholes_gamma():
    mg = mixed_gamma(S, K, T, R, SIG, OptionType.CALL, n_paths=800_000, seed=1)
    assert mg.price == pytest.approx(bs_gamma(S, K, T, R, SIG),
                                     abs=3.0 * mg.std_error + 2e-4)


@pytest.mark.slow
def test_lower_variance_than_lr_gamma():
    mg = mixed_gamma(S, K, T, R, SIG, OptionType.CALL, n_paths=400_000, seed=1)
    lr = lr_greeks(S, K, T, R, SIG, OptionType.CALL, n_paths=400_000, seed=1)
    assert mg.std_error < lr["gamma_se"]


def test_put_gamma_equals_call_gamma():
    c = mixed_gamma(S, K, T, R, SIG, OptionType.CALL, n_paths=200_000, seed=2)
    p = mixed_gamma(S, K, T, R, SIG, OptionType.PUT, n_paths=200_000, seed=2)
    assert c.price == pytest.approx(p.price, abs=3e-3)


def test_gamma_positive():
    mg = mixed_gamma(S, K, T, R, SIG, OptionType.CALL, n_paths=100_000, seed=3)
    assert mg.price > 0.0


def test_bad_params_raise():
    with pytest.raises(ValueError):
        mixed_gamma(-1, K, T, R, SIG)
