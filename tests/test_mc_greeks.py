"""Monte Carlo Greeks by likelihood-ratio and pathwise methods vs Black-Scholes."""

import pytest

from quantforge import (
    OptionType,
    lr_greeks,
    pathwise_delta,
    lr_digital_delta,
    delta as bs_delta,
    gamma as bs_gamma,
    vega as bs_vega,
    cash_or_nothing,
)


S, K, T, R, SIGMA = 100.0, 100.0, 1.0, 0.05, 0.2


@pytest.mark.slow
def test_lr_greeks_match_black_scholes():
    g = lr_greeks(S, K, T, R, SIGMA, OptionType.CALL, n_paths=400_000, seed=1)
    assert g["delta"] == pytest.approx(bs_delta(S, K, T, R, SIGMA),
                                       abs=3.0 * g["delta_se"] + 5e-3)
    assert g["gamma"] == pytest.approx(bs_gamma(S, K, T, R, SIGMA), abs=2e-3)
    assert g["vega"] == pytest.approx(bs_vega(S, K, T, R, SIGMA),
                                      abs=3.0 * g["vega_se"] + 0.5)


@pytest.mark.slow
def test_pathwise_delta_matches_and_beats_lr_variance():
    pw = pathwise_delta(S, K, T, R, SIGMA, OptionType.CALL,
                        n_paths=400_000, seed=2)
    g = lr_greeks(S, K, T, R, SIGMA, OptionType.CALL, n_paths=400_000, seed=2)
    assert pw.price == pytest.approx(bs_delta(S, K, T, R, SIGMA), abs=5e-3)
    # Pathwise has lower variance than LR for this smooth payoff.
    assert pw.std_error < g["delta_se"]


@pytest.mark.slow
def test_lr_digital_delta_matches_analytic():
    h = 0.01
    d_an = (cash_or_nothing(S + h, K, T, R, SIGMA, OptionType.CALL)
            - cash_or_nothing(S - h, K, T, R, SIGMA, OptionType.CALL)) / (2 * h)
    dd = lr_digital_delta(S, K, T, R, SIGMA, OptionType.CALL,
                          n_paths=800_000, seed=3)
    assert dd.price == pytest.approx(d_an, abs=3.0 * dd.std_error + 1e-4)


def test_put_pathwise_delta_negative():
    pw = pathwise_delta(S, K, T, R, SIGMA, OptionType.PUT, n_paths=50_000, seed=4)
    assert pw.price < 0.0


def test_lr_price_matches_call_price():
    from quantforge import call_price
    g = lr_greeks(S, K, T, R, SIGMA, OptionType.CALL, n_paths=100_000, seed=5)
    assert g["price"] == pytest.approx(call_price(S, K, T, R, SIGMA),
                                       abs=3.0 * g["price_se"] + 1e-2)
