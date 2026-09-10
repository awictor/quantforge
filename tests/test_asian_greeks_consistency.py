"""Turnbull-Wakeman arithmetic-Asian Greeks vs an independent Monte Carlo.

The existing Asian-Greeks tests only check that the reported Greeks match a
finite difference of the *same* closed form. Here we cross-check the
Turnbull-Wakeman analytic Greeks against bump Greeks of the arithmetic-Asian
Monte Carlo (common random numbers), an independent pricing method. TW is a
two-moment approximation, so agreement is expected to a small tolerance rather
than exactly.
"""

import pytest

from quantforge import asian_greeks, arithmetic_asian_mc, OptionType


S, K, T, R, SIGMA = 100.0, 100.0, 1.0, 0.05, 0.3
NP, NS, SEED = 200_000, 100, 42


def _mc(s=S, sigma=SIGMA):
    # Common random numbers (same seed) so the bump differences are low-noise.
    return arithmetic_asian_mc(s, K, T, R, sigma, OptionType.CALL,
                               n_steps=NS, n_paths=NP, seed=SEED).price


@pytest.mark.slow
def test_tw_delta_matches_mc_bump():
    g = asian_greeks(S, K, T, R, SIGMA, OptionType.CALL, average="arithmetic")
    h = 0.5
    mc_delta = (_mc(s=S + h) - _mc(s=S - h)) / (2 * h)
    assert g["delta"] == pytest.approx(mc_delta, abs=0.02)


@pytest.mark.slow
def test_tw_gamma_matches_mc_bump():
    g = asian_greeks(S, K, T, R, SIGMA, OptionType.CALL, average="arithmetic")
    h = 0.5
    mc_gamma = (_mc(s=S + h) - 2 * _mc() + _mc(s=S - h)) / (h * h)
    assert g["gamma"] == pytest.approx(mc_gamma, abs=5e-3)


@pytest.mark.slow
def test_tw_vega_matches_mc_bump():
    g = asian_greeks(S, K, T, R, SIGMA, OptionType.CALL, average="arithmetic")
    hv = 1e-3
    mc_vega = (_mc(sigma=SIGMA + hv) - _mc(sigma=SIGMA - hv)) / (2 * hv)
    # Vega is O(20); allow ~2% for the moment-matching approximation.
    assert g["vega"] == pytest.approx(mc_vega, rel=0.02)


@pytest.mark.slow
def test_tw_price_matches_mc():
    g = asian_greeks(S, K, T, R, SIGMA, OptionType.CALL, average="arithmetic")
    assert g["price"] == pytest.approx(_mc(), abs=0.1)


def test_arithmetic_greeks_bracket_geometric_at_low_vol():
    # At low vol the two averages nearly coincide, so their Greeks converge.
    lo = 0.05
    geo = asian_greeks(S, K, T, R, lo, OptionType.CALL, average="geometric")
    ari = asian_greeks(S, K, T, R, lo, OptionType.CALL, average="arithmetic")
    assert ari["delta"] == pytest.approx(geo["delta"], abs=0.01)
    assert ari["vega"] == pytest.approx(geo["vega"], rel=0.05)
