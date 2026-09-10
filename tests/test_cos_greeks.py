"""Analytic COS-method Greeks vs Black-Scholes and finite differences."""

import pytest

from quantforge import (
    OptionType,
    cos_greeks,
    cos_price,
    delta as bs_delta,
    gamma as bs_gamma,
)
from quantforge.cgmy import _cgmy_psi
from quantforge.nig import _nig_psi


def _gbm(sigma):
    return lambda u: -0.5 * sigma * sigma * u * u


def test_gbm_greeks_match_black_scholes():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    g = cos_greeks(S, K, t, r, 0.0, _gbm(sigma), OptionType.CALL, n_terms=256)
    assert g["delta"] == pytest.approx(bs_delta(S, K, t, r, sigma), abs=1e-6)
    assert g["gamma"] == pytest.approx(bs_gamma(S, K, t, r, sigma), abs=1e-6)


def test_gbm_put_delta_matches_black_scholes():
    S, K, t, r, sigma = 100, 105, 1.0, 0.05, 0.25
    g = cos_greeks(S, K, t, r, 0.0, _gbm(sigma), OptionType.PUT, n_terms=256)
    assert g["delta"] == pytest.approx(
        bs_delta(S, K, t, r, sigma, OptionType.PUT), abs=1e-6)
    # Put and call gamma are equal.
    assert g["gamma"] == pytest.approx(bs_gamma(S, K, t, r, sigma), abs=1e-6)


@pytest.mark.parametrize("K", [90, 100, 110])
def test_cgmy_greeks_match_finite_difference(K):
    S, t, r = 100, 1.0, 0.05
    psi = lambda u: _cgmy_psi(u, 4.0, 5.0, 10.0, 0.5)
    g = cos_greeks(S, K, t, r, 0.0, psi, OptionType.CALL, n_terms=512, L=12)
    h = 0.01

    def px(s):
        return cos_price(s, K, t, r, 0.0, psi, OptionType.CALL, n_terms=512, L=12)

    d_fd = (px(S + h) - px(S - h)) / (2 * h)
    g_fd = (px(S + h) - 2 * px(S) + px(S - h)) / (h * h)
    assert g["delta"] == pytest.approx(d_fd, abs=1e-4)
    assert g["gamma"] == pytest.approx(g_fd, abs=1e-4)


def test_nig_greeks_match_finite_difference():
    S, K, t, r = 100, 100, 1.0, 0.03
    psi = lambda u: _nig_psi(u, 15.0, -5.0, 0.5)
    g = cos_greeks(S, K, t, r, 0.0, psi, OptionType.CALL, n_terms=512, L=12)
    h = 0.01

    def px(s):
        return cos_price(s, K, t, r, 0.0, psi, OptionType.CALL, n_terms=512, L=12)

    assert g["delta"] == pytest.approx((px(S + h) - px(S - h)) / (2 * h), abs=1e-4)


def test_price_field_matches_cos_price():
    S, K, t, r = 100, 100, 1.0, 0.05
    psi = _gbm(0.2)
    g = cos_greeks(S, K, t, r, 0.0, psi, OptionType.CALL, n_terms=256)
    assert g["price"] == pytest.approx(
        cos_price(S, K, t, r, 0.0, psi, OptionType.CALL, n_terms=256), abs=1e-10)


def test_call_delta_in_unit_interval():
    g = cos_greeks(100, 100, 1.0, 0.05, 0.0, _gbm(0.2), OptionType.CALL)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0


def test_requires_positive_time():
    with pytest.raises(ValueError):
        cos_greeks(100, 100, 0.0, 0.05, 0.0, _gbm(0.2))
