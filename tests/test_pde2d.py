"""Two-asset ADI (Peaceman-Rachford) PDE vs Margrabe and Kirk."""

import pytest

from quantforge import (
    adi_spread_option,
    adi_two_asset,
    exchange_option,
    spread_option,
)


S1, S2, T = 100.0, 100.0, 1.0
SIG1, SIG2, RHO = 0.2, 0.25, 0.3
R = 0.05


@pytest.mark.slow
def test_zero_strike_spread_matches_margrabe():
    adi = adi_spread_option(S1, S2, 0.0, T, R, SIG1, SIG2, RHO,
                            n1=80, n2=80, n_time=60)
    mar = exchange_option(S1, S2, T, SIG1, SIG2, RHO)
    assert adi == pytest.approx(mar, abs=1e-2)


@pytest.mark.slow
def test_spread_matches_kirk_approximation():
    adi = adi_spread_option(100, 95, 5.0, T, R, SIG1, SIG2, RHO,
                            n1=80, n2=80, n_time=60)
    kirk = spread_option(100, 95, 5.0, T, R, SIG1, SIG2, RHO)
    assert adi == pytest.approx(kirk, abs=2e-2)


def test_higher_correlation_lowers_spread_price():
    vs = [adi_spread_option(100, 100, 10.0, T, R, 0.2, 0.25, rho,
                            n1=50, n2=50, n_time=40)
          for rho in (-0.5, 0.0, 0.9)]
    assert vs[0] > vs[1] > vs[2]


def test_adi_two_asset_generic_payoff_positive():
    # A max-of-two call payoff should be positive and finite.
    v = adi_two_asset(lambda s1, s2: max(max(s1, s2) - 100, 0.0),
                      S1, S2, T, R, SIG1, SIG2, RHO, n1=50, n2=50, n_time=40)
    import math
    assert v > 0.0 and math.isfinite(v)


def test_zero_time_is_payoff():
    v = adi_spread_option(110, 100, 5.0, 0.0, R, SIG1, SIG2, RHO)
    assert v == pytest.approx(max(110 - 100 - 5, 0.0))


def test_bad_prices_raise():
    with pytest.raises(ValueError):
        adi_spread_option(-1, 100, 5.0, T, R, SIG1, SIG2, RHO)
