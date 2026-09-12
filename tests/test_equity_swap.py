"""Equity swaps: total-return and dividend swaps."""

import pytest

from quantforge import (
    total_return_leg, financing_leg, total_return_swap_value, trs_fair_spread,
    dividend_swap_fair_strike, dividend_swap_value,
    variance_swap_payoff, vega_notional_to_variance_notional,
    volatility_swap_payoff, variance_swap_mtm,
)


DIVS = [2.0, 2.1, 2.2]
DFS = [0.97, 0.94, 0.91]


def test_total_return_leg():
    assert total_return_leg(1e6, 100, 110, 3) == pytest.approx(1e6 * (13 / 100))


def test_financing_leg():
    assert financing_leg(1e6, 0.03, 0.005, 1.0) == pytest.approx(1e6 * 0.035)


def test_trs_value_is_equity_minus_financing():
    assert total_return_swap_value(1e6, 100, 110, 3, 0.03, 0.005, 1.0) == \
        pytest.approx(1e6 * 0.13 - 1e6 * 0.035)


def test_fair_spread_zeroes_swap():
    sp = trs_fair_spread(100, 110, 3, 0.03, 1.0)
    assert total_return_swap_value(1e6, 100, 110, 3, 0.03, sp, 1.0) == \
        pytest.approx(0.0, abs=1e-6)


def test_dividend_fair_strike_is_pv():
    K = dividend_swap_fair_strike(DIVS, DFS)
    assert K == pytest.approx(sum(DIVS[i] * DFS[i] for i in range(3)))


def test_dividend_swap_zero_at_fair_strike():
    K = dividend_swap_fair_strike(DIVS, DFS)
    assert dividend_swap_value(DIVS, K, DFS) == pytest.approx(0.0, abs=1e-9)


def test_dividend_swap_gains_when_realized_exceeds():
    K = dividend_swap_fair_strike(DIVS, DFS)
    assert dividend_swap_value([3, 3, 3], K, DFS) > 0


def test_variance_swap_payoff():
    assert variance_swap_payoff(0.25, 0.20, 1e6) == pytest.approx(1e6 * (0.0625 - 0.04))
    assert variance_swap_payoff(0.20, 0.20, 1e6) == pytest.approx(0.0)


def test_vega_to_variance_notional():
    assert vega_notional_to_variance_notional(50000, 0.20) == pytest.approx(
        50000 / (2 * 0.20))


def test_volatility_swap_payoff_linear():
    assert volatility_swap_payoff(0.25, 0.20, 50000) == pytest.approx(50000 * 0.05)


def test_variance_swap_convexity_dominates_vol_swap():
    K, vega = 0.20, 50000
    vn = vega_notional_to_variance_notional(vega, K)
    # Variance swap pays more than the vol swap both above and below the strike.
    assert variance_swap_payoff(0.30, K, vn) > volatility_swap_payoff(0.30, K, vega)
    assert variance_swap_payoff(0.10, K, vn) > volatility_swap_payoff(0.10, K, vega)


def test_variance_swap_mtm_endpoints():
    incep = variance_swap_mtm(0.0, 0.0625, 0.0, 1.0, 0.20, 1e6, 0.95)
    assert incep == pytest.approx(0.95 * 1e6 * (0.0625 - 0.04))
    expiry = variance_swap_mtm(0.09, 0.0, 1.0, 1.0, 0.20, 1e6, 1.0)
    assert expiry == pytest.approx(1e6 * (0.09 - 0.04))


def test_variance_swap_validation():
    with pytest.raises(ValueError):
        vega_notional_to_variance_notional(50000, 0)
    with pytest.raises(ValueError):
        variance_swap_mtm(0.04, 0.04, 2.0, 1.0, 0.2, 1e6, 1.0)


def test_validation():
    with pytest.raises(ValueError):
        total_return_leg(1e6, 0, 110, 3)
    with pytest.raises(ValueError):
        dividend_swap_fair_strike([2], [0.9, 0.8])
    with pytest.raises(ValueError):
        trs_fair_spread(100, 110, 3, 0.03, 0)
