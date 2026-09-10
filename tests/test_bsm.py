"""Correctness gate for QuantForge.

These tests pin the engine against published reference values, enforce
analytic identities (put-call parity, Greek relationships), and verify every
Greek against a central finite-difference of the price function.
"""

import math

import pytest

from quantforge import (
    price, call_price, put_price,
    delta, gamma, vega, theta, rho, greeks,
    implied_volatility, american_price, OptionType,
)
from quantforge.mathfns import norm_cdf, norm_ppf


# --- Reference values (Hull, Options Futures and Other Derivatives) ---
# S=42, K=40, r=0.10, sigma=0.20, t=0.5 -> call 4.7594, put 0.8086
def test_hull_reference_call():
    c = call_price(42, 40, 0.5, 0.10, 0.20)
    assert c == pytest.approx(4.7594, abs=1e-4)


def test_hull_reference_put():
    p = put_price(42, 40, 0.5, 0.10, 0.20)
    assert p == pytest.approx(0.8086, abs=1e-4)


# --- Normal distribution function accuracy ---
def test_norm_cdf_known_points():
    assert norm_cdf(0.0) == pytest.approx(0.5, abs=1e-15)
    assert norm_cdf(1.96) == pytest.approx(0.9750021049, abs=1e-9)
    assert norm_cdf(-1.96) == pytest.approx(0.0249978951, abs=1e-9)


def test_norm_ppf_inverts_cdf():
    for p in (0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99):
        assert norm_cdf(norm_ppf(p)) == pytest.approx(p, abs=1e-12)


# --- Put-call parity: C - P = S*e^{(b-r)t} - K*e^{-rt} ---
@pytest.mark.parametrize("S,K,t,r,sigma,b", [
    (100, 100, 1.0, 0.05, 0.2, 0.05),
    (90, 100, 0.5, 0.03, 0.35, 0.01),
    (120, 100, 2.0, 0.08, 0.15, 0.0),
    (50, 55, 0.25, 0.02, 0.5, -0.01),
])
def test_put_call_parity(S, K, t, r, sigma, b):
    c = price(S, K, t, r, sigma, OptionType.CALL, b=b)
    p = price(S, K, t, r, sigma, OptionType.PUT, b=b)
    lhs = c - p
    rhs = S * math.exp((b - r) * t) - K * math.exp(-r * t)
    assert lhs == pytest.approx(rhs, abs=1e-10)


# --- Greeks vs central finite differences ---
BASE = dict(S=100.0, K=105.0, t=0.75, r=0.04, sigma=0.25)


@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
def test_delta_fd(ot):
    h = 1e-4
    up = price(BASE["S"] + h, BASE["K"], BASE["t"], BASE["r"], BASE["sigma"], ot)
    dn = price(BASE["S"] - h, BASE["K"], BASE["t"], BASE["r"], BASE["sigma"], ot)
    fd = (up - dn) / (2 * h)
    assert delta(**BASE, option_type=ot) == pytest.approx(fd, abs=1e-6)


@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
def test_gamma_fd(ot):
    h = 1e-3
    up = price(BASE["S"] + h, BASE["K"], BASE["t"], BASE["r"], BASE["sigma"], ot)
    mid = price(BASE["S"], BASE["K"], BASE["t"], BASE["r"], BASE["sigma"], ot)
    dn = price(BASE["S"] - h, BASE["K"], BASE["t"], BASE["r"], BASE["sigma"], ot)
    fd = (up - 2 * mid + dn) / (h * h)
    assert gamma(**BASE) == pytest.approx(fd, abs=1e-4)


@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
def test_vega_fd(ot):
    h = 1e-5
    up = price(BASE["S"], BASE["K"], BASE["t"], BASE["r"], BASE["sigma"] + h, ot)
    dn = price(BASE["S"], BASE["K"], BASE["t"], BASE["r"], BASE["sigma"] - h, ot)
    fd = (up - dn) / (2 * h)
    assert vega(**BASE) == pytest.approx(fd, abs=1e-4)


@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
def test_theta_fd(ot):
    h = 1e-5
    # Calendar theta = -dPrice/dt_expiry. FD bumps time-to-expiry, so negate.
    up = price(BASE["S"], BASE["K"], BASE["t"] + h, BASE["r"], BASE["sigma"], ot)
    dn = price(BASE["S"], BASE["K"], BASE["t"] - h, BASE["r"], BASE["sigma"], ot)
    fd = -(up - dn) / (2 * h)
    assert theta(**BASE, option_type=ot) == pytest.approx(fd, abs=1e-3)


@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
def test_rho_fd(ot):
    h = 1e-6
    # Plain stock case b=r: bump both together.
    up = price(BASE["S"], BASE["K"], BASE["t"], BASE["r"] + h, BASE["sigma"], ot,
               b=BASE["r"] + h)
    dn = price(BASE["S"], BASE["K"], BASE["t"], BASE["r"] - h, BASE["sigma"], ot,
               b=BASE["r"] - h)
    fd = (up - dn) / (2 * h)
    assert rho(**BASE, option_type=ot) == pytest.approx(fd, abs=1e-3)


# --- Implied vol round-trip ---
@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
@pytest.mark.parametrize("sigma", [0.05, 0.15, 0.4, 0.9])
@pytest.mark.parametrize("K", [80, 100, 130])
def test_implied_vol_roundtrip(ot, sigma, K):
    S, t, r = 100.0, 1.0, 0.03
    p = price(S, K, t, r, sigma, ot)
    iv = implied_volatility(p, S, K, t, r, ot)
    # Deep-OTM low-vol options carry almost no vega, so vol is ill-conditioned
    # there; assert the recovered *price* matches instead of the raw vol.
    assert price(S, K, t, r, iv, ot) == pytest.approx(p, abs=1e-8)


def test_implied_vol_rejects_arbitrage():
    with pytest.raises(ValueError):
        # A call cannot be worth more than the (carry-adjusted) spot.
        implied_volatility(200.0, 100, 100, 1.0, 0.05, OptionType.CALL)


# --- American vs European cross-checks ---
def test_american_call_no_dividend_equals_european():
    # Without dividends (b=r), American call == European call.
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    eu = call_price(S, K, t, r, sigma)
    am = american_price(S, K, t, r, sigma, OptionType.CALL, steps=800)
    assert am == pytest.approx(eu, abs=1e-2)


def test_american_put_premium_is_positive():
    # American put should be worth at least its European counterpart.
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    eu = put_price(S, K, t, r, sigma)
    am = american_price(S, K, t, r, sigma, OptionType.PUT, steps=800)
    assert am >= eu - 1e-6
    assert am > eu  # early exercise has value here


def test_greeks_bundle_matches_individual():
    g = greeks(**BASE, option_type=OptionType.CALL)
    assert g.price == pytest.approx(price(**BASE, option_type=OptionType.CALL))
    assert g.delta == pytest.approx(delta(**BASE, option_type=OptionType.CALL))
    assert g.gamma == pytest.approx(gamma(**BASE))
    assert g.vega == pytest.approx(vega(**BASE))


# --- Degenerate limits ---
def test_zero_time_is_intrinsic():
    assert call_price(110, 100, 0.0, 0.05, 0.2) == pytest.approx(10.0)
    assert put_price(90, 100, 0.0, 0.05, 0.2) == pytest.approx(10.0)


def test_zero_vol_is_discounted_intrinsic():
    # Forward = S*e^{rt}, call payoff discounted.
    S, K, t, r = 100, 90, 1.0, 0.05
    fwd = S * math.exp(r * t)
    expected = math.exp(-r * t) * max(fwd - K, 0.0)
    assert call_price(S, K, t, r, 0.0) == pytest.approx(expected, abs=1e-10)
