"""Correlation implied by a Bjerksund-Stensland spread price (implied_spread_correlation_bs)."""

import pytest

from quantforge import (
    spread_option_bs, implied_spread_correlation_bs, OptionType,
)


T, R = 1.0, 0.03
S1, S2, K, SIG1, SIG2 = 100.0, 95.0, 5.0, 0.3, 0.35


@pytest.mark.parametrize("rho", [-0.8, -0.3, 0.0, 0.4, 0.85])
def test_round_trip_recovers_rho(rho):
    px = spread_option_bs(S1, S2, K, T, R, SIG1, SIG2, rho)
    got = implied_spread_correlation_bs(px, S1, S2, K, T, R, SIG1, SIG2)
    assert got == pytest.approx(rho, abs=1e-6)


def test_round_trip_put():
    rho = 0.25
    px = spread_option_bs(S1, S2, K, T, R, SIG1, SIG2, rho,
                          option_type=OptionType.PUT)
    got = implied_spread_correlation_bs(px, S1, S2, K, T, R, SIG1, SIG2,
                                        option_type=OptionType.PUT)
    assert got == pytest.approx(rho, abs=1e-6)


def test_price_outside_range_raises():
    # Above the rho=-1 price (the maximum) there is no implied correlation.
    p_max = spread_option_bs(S1, S2, K, T, R, SIG1, SIG2, -0.999999)
    with pytest.raises(ValueError):
        implied_spread_correlation_bs(p_max + 1.0, S1, S2, K, T, R, SIG1, SIG2)


def test_monotone_price_decreasing_in_rho():
    lo = spread_option_bs(S1, S2, K, T, R, SIG1, SIG2, -0.5)
    hi = spread_option_bs(S1, S2, K, T, R, SIG1, SIG2, 0.5)
    assert lo > hi
