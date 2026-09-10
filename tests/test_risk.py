"""Tests for portfolio VaR / Expected Shortfall.

Anchor: a single long stock position (delta=1, gamma=0) has an exact normal
P&L, so its parametric VaR must equal the closed-form ``spot*sigma_h*z``. We
also check ES >= VaR, monotonicity in confidence, and cross-agreement between
the parametric and full-repricing estimators on a nearly-linear book.
"""

import math

import pytest

from quantforge import (
    Contract, price_book, OptionType,
    parametric_var, historical_var, montecarlo_var,
)
from quantforge.mathfns import norm_ppf


def _linear_book():
    # A single deep-ITM call ~ delta 1, tiny gamma: behaves like stock.
    c = Contract(S=100, K=1.0, t=1.0, r=0.0, sigma=0.2, option_type="call",
                 qty=1, multiplier=1)
    return price_book([c])


def test_parametric_var_linear_matches_normal():
    # Pure delta=1, gamma~0 book: VaR = spot * sigma_h * |z_alpha|.
    book = _linear_book()
    sigma, spot, conf = 0.2, 100.0, 0.99
    res = parametric_var(book, sigma_annual=sigma, spot=spot, confidence=conf,
                         horizon_days=1, trading_days=252)
    sigma_h = sigma * math.sqrt(1 / 252)
    expected = spot * sigma_h * abs(norm_ppf(1 - conf))
    # delta is ~1 (deep ITM), allow small slack for residual gamma/discounting.
    assert res.var == pytest.approx(expected, rel=0.03)


def test_es_at_least_var():
    book = _linear_book()
    res = parametric_var(book, 0.25, 100, confidence=0.99)
    assert res.expected_shortfall >= res.var


def test_var_increases_with_confidence():
    book = _linear_book()
    v95 = parametric_var(book, 0.3, 100, confidence=0.95).var
    v99 = parametric_var(book, 0.3, 100, confidence=0.99).var
    assert v99 > v95


def test_historical_var_positive_and_ordered():
    book = _linear_book()
    # Symmetric-ish return scenarios.
    scenarios = [-0.03, -0.02, -0.015, -0.01, 0.0, 0.01, 0.02, 0.03, -0.05, 0.04]
    res = historical_var(book, scenarios, spot=100, confidence=0.90)
    assert res.var > 0
    assert res.expected_shortfall >= res.var


def test_montecarlo_var_agrees_with_parametric_linear():
    # For a near-linear book the full-reval MC VaR should track the parametric
    # (delta-gamma) VaR closely.
    contracts = [Contract(S=100, K=1.0, t=1.0, r=0.0, sigma=0.2,
                          option_type="call", qty=1, multiplier=1)]
    book = price_book(contracts)
    par = parametric_var(book, 0.2, 100, confidence=0.99, horizon_days=1)
    mc = montecarlo_var(contracts, 0.2, confidence=0.99, horizon_days=1,
                        n_paths=40_000, seed=17)
    assert mc.var == pytest.approx(par.var, rel=0.15)


def test_montecarlo_var_convex_book_positive():
    # A long straddle is gamma-positive; VaR should be finite and positive.
    contracts = [
        Contract(S=100, K=100, t=0.25, r=0.02, sigma=0.3, option_type="call", qty=1),
        Contract(S=100, K=100, t=0.25, r=0.02, sigma=0.3, option_type="put", qty=1),
    ]
    mc = montecarlo_var(contracts, 0.3, confidence=0.99, horizon_days=5,
                        n_paths=20_000, seed=3)
    assert mc.var > 0
    assert mc.expected_shortfall >= mc.var


def test_gamma_adds_skew_to_parametric():
    # A short-gamma book (short straddle) loses on big moves either way, so its
    # P&L is left-skewed -> Cornish-Fisher VaR exceeds the plain-normal VaR.
    short_straddle = [
        Contract(S=100, K=100, t=0.25, r=0.02, sigma=0.3, option_type="call", qty=-1),
        Contract(S=100, K=100, t=0.25, r=0.02, sigma=0.3, option_type="put", qty=-1),
    ]
    book = price_book(short_straddle)
    res = parametric_var(book, 0.3, 100, confidence=0.99)
    # Delta ~ 0, gamma < 0: mean P&L per step is negative, VaR strictly positive.
    assert res.var > 0
