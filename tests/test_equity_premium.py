"""Martin (2013) equity-premium lower bound from the simple-variance index."""

import math

import pytest

from quantforge import equity_premium_lower_bound as epb, svix_from_smile


S0, T, R = 100.0, 1.0, 0.03


def test_flat_vol_bound_near_rf_sigma_squared():
    b = epb(S0, T, R, lambda K: 0.2, n_strikes=301, width=8.0)
    assert b == pytest.approx(math.exp(R * T) * 0.04, abs=2e-3)


def test_higher_vol_raises_bound():
    lo = epb(S0, T, R, lambda K: 0.2, n_strikes=301, width=8.0)
    hi = epb(S0, T, R, lambda K: 0.4, n_strikes=301, width=8.0)
    assert hi > lo


def test_bound_equals_rf_times_svix_variance():
    var, _svix = svix_from_smile(S0, T, R, lambda K: 0.2,
                                 n_strikes=301, width=8.0)
    b = epb(S0, T, R, lambda K: 0.2, n_strikes=301, width=8.0)
    assert b == pytest.approx(math.exp(R * T) * var, abs=1e-12)


def test_bound_positive():
    assert epb(S0, T, R, lambda K: 0.25, n_strikes=201) > 0.0


def test_bad_tenor_raises():
    with pytest.raises(ValueError):
        epb(S0, 0.0, R, lambda K: 0.2)
