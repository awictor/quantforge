"""Implied correlation from a Margrabe exchange-option price."""

import pytest

from quantforge import exchange_option, implied_exchange_correlation


S1, S2, T, SIG1, SIG2 = 100.0, 95.0, 1.0, 0.2, 0.25


@pytest.mark.parametrize("rho", [-0.5, 0.0, 0.3, 0.7])
def test_round_trip(rho):
    price = exchange_option(S1, S2, T, SIG1, SIG2, rho)
    assert implied_exchange_correlation(price, S1, S2, T, SIG1, SIG2) == pytest.approx(
        rho, abs=1e-6)


def test_price_monotone_decreasing_in_rho():
    lo = exchange_option(S1, S2, T, SIG1, SIG2, -0.5)
    mid = exchange_option(S1, S2, T, SIG1, SIG2, 0.0)
    hi = exchange_option(S1, S2, T, SIG1, SIG2, 0.5)
    assert lo > mid > hi


def test_with_dividends_round_trip():
    price = exchange_option(S1, S2, T, SIG1, SIG2, 0.4, q1=0.02, q2=0.01)
    rec = implied_exchange_correlation(price, S1, S2, T, SIG1, SIG2, q1=0.02, q2=0.01)
    assert rec == pytest.approx(0.4, abs=1e-6)


def test_out_of_range_raises():
    with pytest.raises(ValueError):
        implied_exchange_correlation(1e-6, S1, S2, T, SIG1, SIG2)
