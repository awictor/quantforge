"""Forward price and put-call parity residual (bsm module)."""

import math

import pytest

from quantforge import (
    forward_price, put_call_parity_residual, call_price, put_price,
)


S, K, T, R, SIG = 100.0, 105.0, 1.0, 0.05, 0.2


def test_forward_price():
    assert forward_price(S, T, R) == pytest.approx(S * math.exp(R * T), abs=1e-12)


def test_forward_future_is_spot():
    assert forward_price(S, T, R, b=0.0) == pytest.approx(S, abs=1e-12)


def test_forward_with_dividend():
    b = R - 0.02
    assert forward_price(S, T, R, b=b) == pytest.approx(S * math.exp(b * T), abs=1e-12)


def test_bsm_prices_satisfy_parity():
    c = call_price(S, K, T, R, SIG)
    p = put_price(S, K, T, R, SIG)
    assert put_call_parity_residual(c, p, S, K, T, R) == pytest.approx(0.0, abs=1e-10)


def test_parity_with_dividend():
    b = R - 0.02
    c = call_price(S, K, T, R, SIG, b=b)
    p = put_price(S, K, T, R, SIG, b=b)
    assert put_call_parity_residual(c, p, S, K, T, R, b=b) == pytest.approx(0.0, abs=1e-10)


def test_mispricing_shows_residual():
    c = call_price(S, K, T, R, SIG)
    p = put_price(S, K, T, R, SIG)
    assert put_call_parity_residual(c + 1.0, p, S, K, T, R) == pytest.approx(1.0, abs=1e-9)
