"""Tests for the Bjerksund-Stensland 1993 single-boundary American pricer."""

import pytest

from quantforge import (
    bjerksund_stensland_1993 as bs93, bjerksund_stensland as bs02,
    american_price, call_price, put_price, OptionType,
)


def test_no_dividend_call_equals_european():
    v = bs93(100, 100, 1.0, 0.05, 0.2, OptionType.CALL)
    assert v == pytest.approx(call_price(100, 100, 1.0, 0.05, 0.2), abs=1e-9)


@pytest.mark.slow
@pytest.mark.parametrize("S,K,t,r,sigma,ot,b", [
    (42, 40, 0.75, 0.04, 0.35, OptionType.CALL, -0.04),
    (100, 100, 0.5, 0.05, 0.30, OptionType.PUT, 0.05),
    (90, 100, 1.0, 0.08, 0.20, OptionType.PUT, 0.08),
])
def test_close_to_binomial_tree(S, K, t, r, sigma, ot, b):
    v = bs93(S, K, t, r, sigma, ot, b=b)
    tree = american_price(S, K, t, r, sigma, ot, b=b, steps=2000)
    # BS1993 is a single-boundary approximation: within ~0.15 of the tree.
    assert v == pytest.approx(tree, abs=0.15)


def test_close_to_bs2002():
    v93 = bs93(100, 100, 0.5, 0.05, 0.3, OptionType.PUT)
    v02 = bs02(100, 100, 0.5, 0.05, 0.3, OptionType.PUT)
    assert v93 == pytest.approx(v02, abs=0.05)


def test_american_put_at_least_european():
    v = bs93(90, 100, 1.0, 0.05, 0.3, OptionType.PUT)
    assert v >= put_price(90, 100, 1.0, 0.05, 0.3) - 1e-6


def test_deep_itm_call_with_dividends_near_intrinsic():
    v = bs93(150, 100, 1.0, 0.05, 0.2, OptionType.CALL, b=-0.10)
    assert v >= (150 - 100) - 1e-6


def test_zero_time_is_intrinsic():
    assert bs93(120, 100, 0.0, 0.05, 0.2, OptionType.CALL) == pytest.approx(20.0)
    assert bs93(80, 100, 0.0, 0.05, 0.2, OptionType.PUT) == pytest.approx(20.0)
