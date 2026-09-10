"""Exact Stulz (1982) rainbow (best-of/worst-of) closed forms."""

import pytest

from quantforge import (
    best_of_call_closed,
    worst_of_call_closed,
    best_of_call,
    worst_of_call,
    call_price,
)


S1, S2, K, T, R = 100.0, 95.0, 100.0, 1.0, 0.05
SIG1, SIG2, RHO = 0.2, 0.3, 0.4


def test_stulz_identity_best_plus_worst():
    # C_max + C_min = c(S1) + c(S2) at the same strike (Stulz).
    bc = best_of_call_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    wc = worst_of_call_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    c1 = call_price(S1, K, T, R, SIG1)
    c2 = call_price(S2, K, T, R, SIG2)
    assert bc + wc == pytest.approx(c1 + c2, abs=1e-9)


@pytest.mark.slow
def test_best_closed_matches_mc():
    bc = best_of_call_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    mc = best_of_call(S1, S2, K, T, R, SIG1, SIG2, RHO, n_paths=600_000, seed=1)
    assert bc == pytest.approx(mc, abs=0.03)


@pytest.mark.slow
def test_worst_closed_matches_mc():
    wc = worst_of_call_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    mc = worst_of_call(S1, S2, K, T, R, SIG1, SIG2, RHO, n_paths=600_000, seed=2)
    assert wc == pytest.approx(mc, abs=0.03)


@pytest.mark.slow
def test_negative_correlation_with_dividends():
    bc = best_of_call_closed(S1, S2, K, T, R, SIG1, SIG2, -0.5, q1=0.02, q2=0.01)
    mc = best_of_call(S1, S2, K, T, R, SIG1, SIG2, -0.5, q1=0.02, q2=0.01,
                      n_paths=1_000_000, seed=3)
    assert bc == pytest.approx(mc, abs=0.05)


def test_best_at_least_worst():
    bc = best_of_call_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    wc = worst_of_call_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    assert bc > wc


def test_best_at_least_max_vanilla():
    # A call on the max is worth at least a vanilla call on either single asset.
    bc = best_of_call_closed(S1, S2, K, T, R, SIG1, SIG2, RHO)
    assert bc >= call_price(S1, K, T, R, SIG1) - 1e-9
    assert bc >= call_price(S2, K, T, R, SIG2) - 1e-9


def test_bad_rho_raises():
    with pytest.raises(ValueError):
        best_of_call_closed(S1, S2, K, T, R, SIG1, SIG2, 1.5)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        worst_of_call_closed(-1, S2, K, T, R, SIG1, SIG2, RHO)
