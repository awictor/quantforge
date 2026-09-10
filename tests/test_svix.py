"""Martin (2013) simple-variance index (SVIX) from a smile."""

import math

import pytest

from quantforge import svix_from_smile, vix_from_smile


S0, T, R = 100.0, 30.0 / 365.0, 0.03


@pytest.mark.parametrize("sig", [0.15, 0.2, 0.4])
def test_flat_smile_near_100_sigma(sig):
    _var, svix = svix_from_smile(S0, T, R, lambda K: sig,
                                 n_strikes=301, width=8.0)
    assert svix == pytest.approx(100.0 * sig, abs=0.2)


def test_svix_is_sqrt_variance():
    var, svix = svix_from_smile(S0, T, R, lambda K: 0.25, n_strikes=201)
    assert svix == pytest.approx(100.0 * math.sqrt(var))


def test_positive():
    _var, svix = svix_from_smile(S0, T, R, lambda K: 0.2, n_strikes=201)
    assert svix > 0.0


def test_differs_from_vix_under_skew():
    # SVIX (1/F^2 weights) and VIX (1/K^2 weights) coincide only to leading
    # order; a skew separates them.
    def down(K):
        return max(0.05, 0.2 + 0.15 * math.log(S0 / K))
    _v, svix = svix_from_smile(S0, T, R, down, n_strikes=301, width=8.0)
    _v2, vix = vix_from_smile(S0, T, R, down, n_strikes=301, width=8.0)
    assert abs(svix - vix) > 1e-3


def test_bad_tenor_raises():
    with pytest.raises(ValueError):
        svix_from_smile(S0, 0.0, R, lambda K: 0.2)
