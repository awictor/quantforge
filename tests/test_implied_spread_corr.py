"""Tests for correlation implied from a spread-option quote."""

import pytest

from quantforge import spread_option, implied_spread_correlation, OptionType


@pytest.mark.parametrize("rho", [-0.5, 0.0, 0.3, 0.8])
def test_round_trip(rho):
    p = spread_option(100, 95, 5, 1.0, 0.05, 0.2, 0.25, rho, option_type=OptionType.CALL)
    back = implied_spread_correlation(p, 100, 95, 5, 1.0, 0.05, 0.2, 0.25,
                                      option_type=OptionType.CALL)
    assert back == pytest.approx(rho, abs=1e-5)


def test_higher_price_lower_correlation():
    # Kirk price decreases in rho, so a richer quote implies a lower correlation.
    lo_price = spread_option(100, 95, 5, 1.0, 0.05, 0.2, 0.25, 0.6)
    hi_price = spread_option(100, 95, 5, 1.0, 0.05, 0.2, 0.25, -0.2)
    rho_from_lo = implied_spread_correlation(lo_price, 100, 95, 5, 1.0, 0.05, 0.2, 0.25)
    rho_from_hi = implied_spread_correlation(hi_price, 100, 95, 5, 1.0, 0.05, 0.2, 0.25)
    assert rho_from_hi < rho_from_lo


def test_out_of_range_price_raises():
    with pytest.raises(ValueError):
        implied_spread_correlation(1000.0, 100, 95, 5, 1.0, 0.05, 0.2, 0.25)


def test_put_round_trip():
    p = spread_option(100, 95, 5, 1.0, 0.05, 0.2, 0.25, 0.4, option_type=OptionType.PUT)
    back = implied_spread_correlation(p, 100, 95, 5, 1.0, 0.05, 0.2, 0.25,
                                      option_type=OptionType.PUT)
    assert back == pytest.approx(0.4, abs=1e-5)
