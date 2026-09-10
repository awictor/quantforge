"""Bermudan/American LSM Greeks by common-random-number bumps."""

import pytest

from quantforge import (
    OptionType,
    bermudan_lsm_greeks,
    bermudan_lsm,
    american_price as crr,
)


S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.2


@pytest.mark.slow
def test_delta_matches_binomial():
    g = bermudan_lsm_greeks(S, K, T, R, SIG, OptionType.PUT,
                            n_steps=50, n_paths=200_000, seed=1, h_rel=0.02)
    h = 0.5
    up = crr(S + h, K, T, R, SIG, OptionType.PUT, b=R, steps=2000)
    dn = crr(S - h, K, T, R, SIG, OptionType.PUT, b=R, steps=2000)
    assert g["delta"] == pytest.approx((up - dn) / (2 * h), abs=1e-2)


@pytest.mark.slow
def test_gamma_positive_and_roughly_right():
    # LSM gamma is a noisy second difference over a re-fit regression -- only
    # indicative, so a loose band around the binomial gamma.
    g = bermudan_lsm_greeks(S, K, T, R, SIG, OptionType.PUT,
                            n_steps=50, n_paths=200_000, seed=1, h_rel=0.02)
    assert g["gamma"] > 0.0
    h = 0.5
    up = crr(S + h, K, T, R, SIG, OptionType.PUT, b=R, steps=2000)
    dn = crr(S - h, K, T, R, SIG, OptionType.PUT, b=R, steps=2000)
    base = crr(S, K, T, R, SIG, OptionType.PUT, b=R, steps=2000)
    tree_gamma = (up - 2 * base + dn) / (h * h)
    assert g["gamma"] == pytest.approx(tree_gamma, abs=2e-2)


def test_put_delta_negative():
    g = bermudan_lsm_greeks(S, K, T, R, SIG, OptionType.PUT,
                            n_steps=40, n_paths=40_000, seed=2)
    assert g["delta"] < 0.0


def test_price_field_matches_lsm():
    g = bermudan_lsm_greeks(S, K, T, R, SIG, OptionType.PUT,
                            n_steps=40, n_paths=30_000, seed=3)
    direct = bermudan_lsm(S, K, T, R, SIG, OptionType.PUT, n_steps=40,
                          n_paths=30_000, seed=3)
    assert g["price"] == pytest.approx(direct, abs=1e-9)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        bermudan_lsm_greeks(-1, K, T, R, SIG)
