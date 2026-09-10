"""Tests for the optional NumPy fast path.

The vectorized results must agree with the scalar engine to machine precision,
and scalar/array broadcasting must work. Skipped entirely if NumPy is absent.
"""

import pytest

np = pytest.importorskip("numpy")

from quantforge import (
    price_array, delta_array, gamma_array, vega_array, greeks_array,
    call_price, put_price, delta, gamma, vega, OptionType,
)


def test_price_array_matches_scalar():
    S = np.array([90.0, 100.0, 110.0])
    K = np.array([100.0, 100.0, 100.0])
    t, r, sigma = 0.5, 0.04, 0.25
    got = price_array(S, K, t, r, sigma, OptionType.CALL)
    for i in range(len(S)):
        exact = call_price(float(S[i]), float(K[i]), t, r, sigma)
        assert got[i] == pytest.approx(exact, abs=1e-12)


def test_put_array_matches_scalar():
    S = np.linspace(80, 120, 25)
    got = price_array(S, 100.0, 1.0, 0.03, 0.3, "put")
    for i, s in enumerate(S):
        assert got[i] == pytest.approx(put_price(float(s), 100.0, 1.0, 0.03, 0.3), abs=1e-12)


def test_scalar_broadcasts_against_strike_grid():
    K = np.linspace(80, 120, 50)
    got = price_array(100.0, K, 0.75, 0.05, 0.2, OptionType.CALL)
    assert got.shape == (50,)
    assert got[0] > got[-1]  # deeper ITM call worth more


def test_greeks_array_matches_scalar():
    S = np.array([95.0, 100.0, 105.0])
    K, t, r, sigma = 100.0, 0.5, 0.04, 0.25
    g = greeks_array(S, K, t, r, sigma, OptionType.CALL)
    for i, s in enumerate(S):
        s = float(s)
        assert g["delta"][i] == pytest.approx(delta(s, K, t, r, sigma, OptionType.CALL), abs=1e-12)
        assert g["gamma"][i] == pytest.approx(gamma(s, K, t, r, sigma), abs=1e-12)
        assert g["vega"][i] == pytest.approx(vega(s, K, t, r, sigma), abs=1e-12)


def test_large_chain_runs():
    S = np.full(10_000, 100.0)
    K = np.linspace(50, 150, 10_000)
    p = price_array(S, K, 1.0, 0.05, 0.3, OptionType.CALL)
    assert p.shape == (10_000,)
    assert np.all(p >= 0)
