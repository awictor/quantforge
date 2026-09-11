"""Implied forward from a call-put pair via parity (bsm module)."""

import pytest

from quantforge import (
    implied_forward_from_parity, forward_price, call_price, put_price,
)


S, K, T, R, SIG = 100.0, 105.0, 1.0, 0.05, 0.2


def test_recovers_forward_no_dividend():
    c = call_price(S, K, T, R, SIG)
    p = put_price(S, K, T, R, SIG)
    assert implied_forward_from_parity(c, p, K, T, R) == pytest.approx(
        forward_price(S, T, R), abs=1e-9)


def test_recovers_forward_with_dividend():
    b = R - 0.03
    c = call_price(S, K, T, R, SIG, b=b)
    p = put_price(S, K, T, R, SIG, b=b)
    assert implied_forward_from_parity(c, p, K, T, R) == pytest.approx(
        forward_price(S, T, R, b=b), abs=1e-9)


@pytest.mark.parametrize("strike", [90.0, 105.0, 120.0])
def test_strike_independent(strike):
    c = call_price(S, strike, T, R, SIG)
    p = put_price(S, strike, T, R, SIG)
    assert implied_forward_from_parity(c, p, strike, T, R) == pytest.approx(
        forward_price(S, T, R), abs=1e-9)


def test_atm_forward_call_equals_put():
    # At the forward strike the call and put are equal, so implied F == K.
    F = forward_price(S, T, R)
    c = call_price(S, F, T, R, SIG)
    p = put_price(S, F, T, R, SIG)
    assert c == pytest.approx(p, abs=1e-9)
    assert implied_forward_from_parity(c, p, F, T, R) == pytest.approx(F, abs=1e-9)
