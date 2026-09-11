"""Implied uniform correlation from a geometric-basket price."""

import pytest

from quantforge import (
    geometric_basket_option, implied_geometric_basket_correlation,
)


S = [100.0, 95.0, 105.0]
W = [0.5, 0.3, 0.2]
SIG = [0.2, 0.25, 0.18]
K, T, R = 100.0, 1.0, 0.05


def _corr(rho, n=3):
    return [[1.0 if i == j else rho for j in range(n)] for i in range(n)]


@pytest.mark.parametrize("rho", [-0.3, 0.0, 0.4, 0.8])
@pytest.mark.parametrize("ot", ["call", "put"])
def test_round_trip(rho, ot):
    target = geometric_basket_option(S, W, K, T, R, SIG, _corr(rho), None, ot)
    rec = implied_geometric_basket_correlation(target, S, W, K, T, R, SIG, None, ot)
    assert rec == pytest.approx(rho, abs=1e-5)


def test_two_asset_round_trip():
    s, w, sig = [100.0, 90.0], [0.6, 0.4], [0.2, 0.3]
    target = geometric_basket_option(s, w, K, T, R, sig, _corr(0.5, 2), None, "call")
    rec = implied_geometric_basket_correlation(target, s, w, K, T, R, sig, None, "call")
    assert rec == pytest.approx(0.5, abs=1e-5)


def test_out_of_range_raises():
    # A price above the rho=1 value cannot be matched.
    hi = geometric_basket_option(S, W, K, T, R, SIG, _corr(0.999999), None, "call")
    with pytest.raises(ValueError):
        implied_geometric_basket_correlation(hi * 1.5, S, W, K, T, R, SIG, None, "call")


def test_needs_two_assets():
    with pytest.raises(ValueError):
        implied_geometric_basket_correlation(5.0, [100.0], [1.0], K, T, R, [0.2],
                                             None, "call")
