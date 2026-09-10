"""Tests for the Bjerksund-Stensland (2002) American approximation.

BS2002 is a closed-form *approximation*, so we validate it against the
converged binomial tree within its known accuracy band (a few cents for
typical inputs) and pin the exact identities it must satisfy: no-dividend
American calls equal European calls, prices dominate intrinsic value, and puts
gain from early exercise versus their European counterpart.
"""

import math

import pytest

from quantforge import (
    bjerksund_stensland as bs, american_price, call_price, put_price,
    OptionType,
)


TREE_CASES = [
    # (S, K, t, r, sigma, type, b)
    (42, 40, 0.75, 0.04, 0.35, OptionType.CALL, -0.04),   # dividend call
    (90, 100, 1.0, 0.08, 0.20, OptionType.PUT, 0.08),
    (100, 100, 0.5, 0.05, 0.30, OptionType.PUT, 0.05),
    (110, 100, 0.5, 0.10, 0.25, OptionType.PUT, 0.10),
    (100, 110, 1.0, 0.03, 0.40, OptionType.PUT, 0.03),
]


@pytest.mark.parametrize("S,K,t,r,sigma,ot,b", TREE_CASES)
def test_bs2002_close_to_binomial_tree(S, K, t, r, sigma, ot, b):
    approx = bs(S, K, t, r, sigma, ot, b=b)
    tree = american_price(S, K, t, r, sigma, ot, b=b, steps=3000)
    # BS2002 is accurate to a few cents against the converged tree.
    assert approx == pytest.approx(tree, abs=0.10)


def test_no_dividend_call_equals_european():
    # b = r (no dividends): early exercise is never optimal, so BS2002 returns
    # the European Black-Scholes value.
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    assert bs(S, K, t, r, sigma, OptionType.CALL) == pytest.approx(
        call_price(S, K, t, r, sigma), abs=1e-9)


def test_american_put_at_least_european():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    am = bs(S, K, t, r, sigma, OptionType.PUT)
    eu = put_price(S, K, t, r, sigma)
    assert am >= eu - 1e-6
    assert am > eu  # early exercise has value


def test_price_dominates_intrinsic():
    # A deep ITM American put is worth at least its immediate-exercise value.
    S, K, t, r, sigma = 60, 100, 0.5, 0.05, 0.3
    am = bs(S, K, t, r, sigma, OptionType.PUT)
    assert am >= (K - S) - 1e-6


def test_deep_itm_call_with_dividends_near_intrinsic():
    # Deep ITM call with heavy dividends -> early exercise -> ~ intrinsic.
    S, K, t, r, sigma = 150, 100, 1.0, 0.05, 0.2
    v = bs(S, K, t, r, sigma, OptionType.CALL, b=-0.10)
    assert v >= (S - K) - 1e-6


def test_zero_time_is_intrinsic():
    assert bs(120, 100, 0.0, 0.05, 0.2, OptionType.CALL) == pytest.approx(20.0)
    assert bs(80, 100, 0.0, 0.05, 0.2, OptionType.PUT) == pytest.approx(20.0)


def test_bivariate_normal_reference_values():
    from quantforge.american import _bivariate_normal as B
    # B(0,0,rho) = 1/4 + asin(rho)/(2 pi).
    for rho in (-0.5, 0.0, 0.5, 0.786):
        assert B(0, 0, rho) == pytest.approx(0.25 + math.asin(rho) / (2 * math.pi), abs=1e-6)
    # Independence: B(a,b,0) = N(a) N(b).
    from quantforge.mathfns import norm_cdf
    assert B(0.5, 0.3, 0.0) == pytest.approx(norm_cdf(0.5) * norm_cdf(0.3), abs=1e-9)
