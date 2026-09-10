"""Kim American-put Greeks (one boundary solve, bumped price) vs binomial."""

import pytest

from quantforge import (
    OptionType,
    kim_put_greeks,
    kim_american_put,
    american_price as crr,
)


CASES = [
    (100, 100, 1.0, 0.05, 0.20, 0.0),
    (95, 100, 0.5, 0.05, 0.30, 0.02),
    (110, 100, 1.0, 0.08, 0.25, 0.04),
]


@pytest.mark.slow
@pytest.mark.parametrize("S,K,t,r,sigma,q", CASES)
def test_delta_matches_binomial_bump(S, K, t, r, sigma, q):
    g = kim_put_greeks(S, K, t, r, sigma, q=q, n_steps=160)
    h = 0.5
    up = crr(S + h, K, t, r, sigma, OptionType.PUT, b=r - q, steps=5000)
    dn = crr(S - h, K, t, r, sigma, OptionType.PUT, b=r - q, steps=5000)
    assert g["delta"] == pytest.approx((up - dn) / (2 * h), abs=5e-3)


@pytest.mark.slow
def test_gamma_matches_binomial_bump():
    S, K, t, r, sigma, q = 100, 100, 1.0, 0.05, 0.2, 0.0
    g = kim_put_greeks(S, K, t, r, sigma, q=q, n_steps=200)
    h = 0.5
    up = crr(S + h, K, t, r, sigma, OptionType.PUT, b=r, steps=6000)
    dn = crr(S - h, K, t, r, sigma, OptionType.PUT, b=r, steps=6000)
    base = crr(S, K, t, r, sigma, OptionType.PUT, b=r, steps=6000)
    tree_gamma = (up - 2 * base + dn) / (h * h)
    assert g["gamma"] == pytest.approx(tree_gamma, abs=2e-3)


def test_put_delta_negative_gamma_positive():
    g = kim_put_greeks(100, 100, 1.0, 0.05, 0.2, n_steps=120)
    assert g["delta"] < 0.0
    assert g["gamma"] > 0.0


def test_price_field_matches_direct():
    g = kim_put_greeks(100, 100, 1.0, 0.05, 0.2, q=0.02, n_steps=120)
    direct = kim_american_put(100, 100, 1.0, 0.05, 0.2, q=0.02, n_steps=120)
    assert g["price"] == pytest.approx(direct, abs=1e-9)


def test_theta_negative_for_put():
    # A vanilla American put loses value as time passes (theta < 0 here).
    g = kim_put_greeks(100, 100, 1.0, 0.05, 0.2, n_steps=120)
    assert g["theta"] < 0.0


def test_bad_params_raise():
    with pytest.raises(ValueError):
        kim_put_greeks(-1, 100, 1.0, 0.05, 0.2)
    with pytest.raises(ValueError):
        kim_put_greeks(100, 100, 0.0, 0.05, 0.2)
