"""Tests for epsilon (dividend rho): dPrice/dq."""

import math

import pytest

from quantforge import epsilon, price, OptionType


CASES = [
    dict(S=100, K=90, t=0.75, r=0.04, sigma=0.25),
    dict(S=100, K=100, t=0.5, r=0.03, sigma=0.3),
    dict(S=100, K=115, t=1.0, r=0.05, sigma=0.2),
]


@pytest.mark.parametrize("p", CASES)
@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
def test_epsilon_matches_finite_difference(p, ot):
    q = 0.03
    b = p["r"] - q
    h = 1e-6
    up = price(p["S"], p["K"], p["t"], p["r"], p["sigma"], ot, b=p["r"] - (q + h))
    dn = price(p["S"], p["K"], p["t"], p["r"], p["sigma"], ot, b=p["r"] - (q - h))
    fd = (up - dn) / (2 * h)
    assert epsilon(**p, option_type=ot, b=b) == pytest.approx(fd, abs=1e-3)


def test_call_epsilon_negative_put_positive():
    # Raising the dividend yield lowers the forward: call loses, put gains.
    b = 0.04 - 0.03
    assert epsilon(100, 100, 1.0, 0.04, 0.25, OptionType.CALL, b=b) < 0
    assert epsilon(100, 100, 1.0, 0.04, 0.25, OptionType.PUT, b=b) > 0


def test_epsilon_zero_at_expiry():
    assert epsilon(100, 100, 0.0, 0.04, 0.25, OptionType.CALL) == 0.0


def test_epsilon_scales_with_time():
    # Longer-dated options have larger dividend sensitivity in magnitude.
    short = abs(epsilon(100, 100, 0.25, 0.04, 0.25, OptionType.CALL, b=0.01))
    long = abs(epsilon(100, 100, 2.0, 0.04, 0.25, OptionType.CALL, b=0.01))
    assert long > short
