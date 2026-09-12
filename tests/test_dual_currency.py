"""Dual-currency deposits."""

import pytest

from quantforge import (
    dcd_enhanced_yield, dcd_option_premium_rate, dcd_maturity_payoff,
    dcd_breakeven_spot,
)


SPOT, STRIKE, T, RD, RF, SIG = 1.10, 1.12, 0.25, 0.05, 0.03, 0.10
BASE = 0.05


def test_premium_positive():
    assert dcd_option_premium_rate(SPOT, STRIKE, T, RD, RF, SIG) > 0


def test_enhanced_yield_above_base():
    pr = dcd_option_premium_rate(SPOT, STRIKE, T, RD, RF, SIG)
    assert dcd_enhanced_yield(BASE, pr, T) > BASE


def test_volatility_raises_premium():
    assert dcd_option_premium_rate(SPOT, STRIKE, T, RD, RF, 0.15) > \
        dcd_option_premium_rate(SPOT, STRIKE, T, RD, RF, 0.10)


def test_not_converted_below_strike():
    cpn = 0.10
    gross = 100000 * (1 + cpn * T)
    assert dcd_maturity_payoff(100000, cpn, T, 1.08, STRIKE) == pytest.approx(gross)


def test_converted_above_strike():
    cpn = 0.10
    gross = 100000 * (1 + cpn * T)
    conv = dcd_maturity_payoff(100000, cpn, T, 1.16, STRIKE)
    assert conv < gross
    assert conv == pytest.approx(gross * STRIKE / 1.16)


def test_breakeven_above_strike_and_matches_plain():
    cpn = 0.10
    bstar = dcd_breakeven_spot(cpn, BASE, T, STRIKE)
    assert bstar > STRIKE
    plain = 100000 * (1 + BASE * T)
    assert dcd_maturity_payoff(100000, cpn, T, bstar, STRIKE) == pytest.approx(plain, abs=1e-6)


def test_validation():
    with pytest.raises(ValueError):
        dcd_option_premium_rate(0, STRIKE, T, RD, RF, SIG)
    with pytest.raises(ValueError):
        dcd_enhanced_yield(BASE, 0.01, 0)
