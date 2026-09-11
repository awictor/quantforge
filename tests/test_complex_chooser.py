"""Complex chooser option (Rubinstein 1991) — chooser.complex_chooser_option."""

import math

import pytest

from quantforge import (
    complex_chooser_option, complex_chooser_option_greeks, chooser_option,
)
from quantforge.bsm import call_price, put_price


S, R, SIG = 100.0, 0.05, 0.25


def test_reduces_to_simple_chooser():
    # Equal strikes and maturities on both legs = the simple chooser.
    tc, T = 0.5, 1.0
    cx = complex_chooser_option(S, 100.0, 100.0, tc, T, T, R, SIG)
    simple = chooser_option(S, 100.0, tc, T, R, SIG)
    assert cx == pytest.approx(simple, abs=1e-8)


def test_matches_monte_carlo():
    Kc, Kp, tc, Tc, Tp = 105.0, 95.0, 0.5, 1.0, 1.5
    cf = complex_chooser_option(S, Kc, Kp, tc, Tc, Tp, R, SIG)
    import random
    rng = random.Random(5)
    N = 200000
    disc = math.exp(-R * tc)
    acc = 0.0
    for _ in range(N):
        z = rng.gauss(0.0, 1.0)
        Stc = S * math.exp((R - 0.5 * SIG * SIG) * tc + SIG * math.sqrt(tc) * z)
        acc += max(call_price(Stc, Kc, Tc - tc, R, SIG),
                   put_price(Stc, Kp, Tp - tc, R, SIG))
    mc = disc * acc / N
    assert cf == pytest.approx(mc, rel=0.01)


def test_more_valuable_than_either_leg_today():
    # The right to choose is worth at least each leg valued outright today.
    Kc, Kp, tc, Tc, Tp = 105.0, 95.0, 0.5, 1.0, 1.5
    cx = complex_chooser_option(S, Kc, Kp, tc, Tc, Tp, R, SIG)
    assert cx > call_price(S, Kc, Tc, R, SIG)
    assert cx > put_price(S, Kp, Tp, R, SIG)


def test_greeks_price_field_and_delta_fd():
    Kc, Kp, tc, Tc, Tp = 105.0, 95.0, 0.5, 1.0, 1.5
    g = complex_chooser_option_greeks(S, Kc, Kp, tc, Tc, Tp, R, SIG)
    assert g["price"] == pytest.approx(
        complex_chooser_option(S, Kc, Kp, tc, Tc, Tp, R, SIG), abs=1e-9)
    h = 1e-4 * S
    fd = (complex_chooser_option(S + h, Kc, Kp, tc, Tc, Tp, R, SIG)
          - complex_chooser_option(S - h, Kc, Kp, tc, Tc, Tp, R, SIG)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-5)
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0


def test_validation():
    with pytest.raises(ValueError):
        complex_chooser_option(S, 105.0, 95.0, 0.0, 1.0, 1.5, R, SIG)  # tc = 0
    with pytest.raises(ValueError):
        complex_chooser_option(S, 105.0, 95.0, 1.2, 1.0, 1.5, R, SIG)  # tc > Tc
