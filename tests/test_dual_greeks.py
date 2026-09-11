"""Dual delta and dual gamma (strike sensitivities) in greeks2."""

import pytest

from quantforge import dual_delta, dual_gamma
from quantforge.bsm import price
from quantforge.exotics import cash_or_nothing


S, K, T, R, SIG = 100.0, 105.0, 1.0, 0.05, 0.2


@pytest.mark.parametrize("ot", ["call", "put"])
def test_dual_delta_matches_fd(ot):
    h = 1e-4
    fd = (price(S, K + h, T, R, SIG, ot) - price(S, K - h, T, R, SIG, ot)) / (2 * h)
    assert dual_delta(S, K, T, R, SIG, ot) == pytest.approx(fd, abs=1e-5)


def test_negative_call_dual_delta_is_cash_digital():
    assert -dual_delta(S, K, T, R, SIG, "call") == pytest.approx(
        cash_or_nothing(S, K, T, R, SIG, "call"), abs=1e-9)


def test_dual_gamma_matches_second_difference():
    h = 1e-4
    fd2 = (price(S, K + h, T, R, SIG, "call") - 2 * price(S, K, T, R, SIG, "call")
           + price(S, K - h, T, R, SIG, "call")) / (h * h)
    assert dual_gamma(S, K, T, R, SIG) == pytest.approx(fd2, abs=1e-3)


def test_dual_gamma_non_negative():
    for k in (80.0, 100.0, 120.0):
        assert dual_gamma(S, k, T, R, SIG) >= 0.0


def test_call_dual_delta_negative_put_positive():
    assert dual_delta(S, K, T, R, SIG, "call") < 0.0
    assert dual_delta(S, K, T, R, SIG, "put") > 0.0


def test_dual_delta_parity():
    # call - put dual delta = -e^{-rt} (N(d2) + N(-d2)) = -e^{-rt}.
    import math
    c = dual_delta(S, K, T, R, SIG, "call")
    p = dual_delta(S, K, T, R, SIG, "put")
    assert c - p == pytest.approx(-math.exp(-R * T), abs=1e-9)
