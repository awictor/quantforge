"""Greeks of the Bjerksund-Stensland spread option (spread_option_bs_greeks)."""

import math

import pytest

from quantforge import (
    spread_option_bs_greeks, spread_option_bs, OptionType,
)


T, R = 1.0, 0.03
S1, S2, K, SIG1, SIG2, RHO = 100.0, 95.0, 5.0, 0.3, 0.3, 0.4


def test_greek_signs():
    g = spread_option_bs_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO)
    assert g["delta1"] > 0.0        # long asset 1
    assert g["delta2"] < 0.0        # short asset 2
    assert g["gamma1"] > 0.0
    assert g["gamma2"] > 0.0
    assert g["cross"] < 0.0         # d2V/dS1dS2 < 0 for a spread
    assert g["corr_vega"] < 0.0     # higher corr -> lower spread vol -> lower call


def test_delta1_matches_independent_bump():
    g = spread_option_bs_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO)
    h = 1e-2
    up = spread_option_bs(S1 + h, S2, K, T, R, SIG1, SIG2, RHO)
    dn = spread_option_bs(S1 - h, S2, K, T, R, SIG1, SIG2, RHO)
    assert g["delta1"] == pytest.approx((up - dn) / (2 * h), abs=1e-4)


def test_put_call_delta_parity():
    # d/dS1 (C - P) = disc * dF1/dS1 = disc * e^{(r)t} = 1 (q1=0); d/dS2 = -1.
    gc = spread_option_bs_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO,
                                 option_type=OptionType.CALL)
    gp = spread_option_bs_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO,
                                 option_type=OptionType.PUT)
    assert gc["delta1"] - gp["delta1"] == pytest.approx(1.0, abs=1e-4)
    assert gc["delta2"] - gp["delta2"] == pytest.approx(-1.0, abs=1e-4)


def test_gamma_matches_cross_magnitude():
    # For a spread of two near-symmetric lognormals the cross-gamma is close to
    # -sqrt(gamma1 * gamma2) (perfect anti-diagonal curvature).
    g = spread_option_bs_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO)
    approx = -math.sqrt(g["gamma1"] * g["gamma2"])
    assert g["cross"] == pytest.approx(approx, rel=0.05)


def test_put_gammas_match_call():
    # Gamma is identical for call and put (parity is linear in spots).
    gc = spread_option_bs_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO,
                                 option_type=OptionType.CALL)
    gp = spread_option_bs_greeks(S1, S2, K, T, R, SIG1, SIG2, RHO,
                                 option_type=OptionType.PUT)
    assert gc["gamma1"] == pytest.approx(gp["gamma1"], abs=1e-6)
    assert gc["gamma2"] == pytest.approx(gp["gamma2"], abs=1e-6)
