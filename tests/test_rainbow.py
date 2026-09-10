"""Tests for best-of / worst-of two-asset (rainbow) options."""

import pytest

from quantforge import best_of_call, worst_of_call, call_price, OptionType


BASE = dict(S1=100, S2=100, K=100, t=1.0, r=0.05, sigma1=0.2, sigma2=0.25, rho=0.4)


@pytest.mark.slow
def test_stulz_identity_best_plus_worst():
    # best-of + worst-of = call(S1) + call(S2) at the same strike (Stulz 1982).
    b = best_of_call(**BASE, n_paths=200_000, seed=1)
    w = worst_of_call(**BASE, n_paths=200_000, seed=1)
    c1 = call_price(100, 100, 1.0, 0.05, 0.2)
    c2 = call_price(100, 100, 1.0, 0.05, 0.25)
    assert b + w == pytest.approx(c1 + c2, abs=0.1)


@pytest.mark.slow
def test_best_dominates_each_single_call():
    b = best_of_call(**BASE, n_paths=100_000, seed=2)
    c1 = call_price(100, 100, 1.0, 0.05, 0.2)
    c2 = call_price(100, 100, 1.0, 0.05, 0.25)
    assert b >= c1 - 0.1
    assert b >= c2 - 0.1


@pytest.mark.slow
def test_worst_below_each_single_call():
    w = worst_of_call(**BASE, n_paths=100_000, seed=3)
    c1 = call_price(100, 100, 1.0, 0.05, 0.2)
    c2 = call_price(100, 100, 1.0, 0.05, 0.25)
    assert w <= c1 + 0.1
    assert w <= c2 + 0.1


def test_best_at_least_worst():
    b = best_of_call(**BASE, n_paths=40_000, seed=4)
    w = worst_of_call(**BASE, n_paths=40_000, seed=4)
    assert b >= w


def test_higher_correlation_widens_gap():
    # As correlation rises the two assets move together, so best and worst
    # converge; low correlation widens the best-minus-worst gap.
    hi = best_of_call(**{**BASE, "rho": 0.95}, n_paths=40_000, seed=5) - \
        worst_of_call(**{**BASE, "rho": 0.95}, n_paths=40_000, seed=5)
    lo = best_of_call(**{**BASE, "rho": -0.5}, n_paths=40_000, seed=5) - \
        worst_of_call(**{**BASE, "rho": -0.5}, n_paths=40_000, seed=5)
    assert lo > hi


def test_reproducible_with_seed():
    a = best_of_call(**BASE, n_paths=10_000, seed=9)
    b = best_of_call(**BASE, n_paths=10_000, seed=9)
    assert a == b
