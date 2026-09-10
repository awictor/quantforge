"""Barone-Adesi-Whaley American approximation vs binomial and Kim."""

import pytest

from quantforge import (
    OptionType,
    baw_american,
    kim_american_put,
    american_price as crr,
    call_price,
    put_price,
)


# (type, S, K, t, r, sigma, q)
CASES = [
    ("put", 100, 100, 1.0, 0.05, 0.20, 0.0),
    ("put", 90, 100, 0.5, 0.05, 0.30, 0.0),
    ("put", 110, 100, 1.0, 0.08, 0.25, 0.04),
    ("call", 100, 100, 1.0, 0.05, 0.20, 0.06),
    ("call", 110, 100, 0.5, 0.03, 0.25, 0.05),
]


@pytest.mark.slow
@pytest.mark.parametrize("typ,S,K,t,r,sigma,q", CASES)
def test_baw_close_to_binomial(typ, S, K, t, r, sigma, q):
    ot = OptionType.CALL if typ == "call" else OptionType.PUT
    baw = baw_american(S, K, t, r, sigma, ot, b=r - q)
    tree = crr(S, K, t, r, sigma, ot, b=r - q, steps=4000)
    # BAW is a quadratic *approximation*: within ~0.1 of the converged tree.
    assert baw == pytest.approx(tree, abs=0.1)


def test_no_dividend_call_equals_european():
    v = baw_american(100, 100, 1.0, 0.05, 0.2, OptionType.CALL, b=0.05)
    assert v == pytest.approx(call_price(100, 100, 1.0, 0.05, 0.2), abs=1e-9)


def test_american_put_at_least_european_and_intrinsic():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    am = baw_american(S, K, t, r, sigma, OptionType.PUT)
    assert am >= put_price(S, K, t, r, sigma) - 1e-6
    assert am >= (K - S) - 1e-6


def test_deep_itm_put_is_intrinsic():
    assert baw_american(50, 100, 1.0, 0.05, 0.2, OptionType.PUT) == pytest.approx(
        50.0, abs=1e-6)


def test_deep_itm_dividend_call_near_intrinsic():
    v = baw_american(160, 100, 1.0, 0.05, 0.2, OptionType.CALL, b=-0.05)
    assert v >= (160 - 100) - 1e-6


@pytest.mark.slow
def test_agrees_with_kim_within_approximation_error():
    # BAW and the (accurate) Kim integral equation should be within BAW's
    # approximation band of each other.
    S, K, t, r, sigma, q = 100, 100, 1.0, 0.05, 0.2, 0.0
    baw = baw_american(S, K, t, r, sigma, OptionType.PUT, b=r - q)
    kim = kim_american_put(S, K, t, r, sigma, q=q, n_steps=200)
    assert baw == pytest.approx(kim, abs=0.05)


def test_zero_time_is_intrinsic():
    assert baw_american(120, 100, 0.0, 0.05, 0.2, OptionType.CALL) == pytest.approx(20.0)
