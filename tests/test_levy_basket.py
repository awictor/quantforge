"""n-asset Levy moment-matched basket option."""

import math

import pytest

from quantforge import levy_basket_option as lb
from quantforge.linalg import basket_option_mc
from quantforge.bsm import call_price


def test_single_asset_equals_bsm():
    p = lb([100], [1.0], 105, 1.0, 0.05, [0.25], [[1.0]])
    assert abs(p - call_price(100, 105, 1.0, 0.05, 0.25)) < 1e-9


def test_three_asset_matches_monte_carlo():
    spots, w, sig = [100, 90, 110], [0.4, 0.3, 0.3], [0.2, 0.25, 0.3]
    corr = [[1, 0.5, 0.3], [0.5, 1, 0.4], [0.3, 0.4, 1]]
    levy = lb(spots, w, 100, 1.0, 0.03, sig, corr)
    m = basket_option_mc(spots, w, 100, 1.0, 0.03, sig, corr, n_paths=200000, seed=1)
    mc = m[0] if isinstance(m, (tuple, list)) else m
    assert abs(levy - mc) < 0.3


def test_put_call_parity():
    spots, w, sig = [100, 90, 110], [0.4, 0.3, 0.3], [0.2, 0.25, 0.3]
    corr = [[1, 0.5, 0.3], [0.5, 1, 0.4], [0.3, 0.4, 1]]
    c = lb(spots, w, 100, 1.0, 0.03, sig, corr, is_call=True)
    p = lb(spots, w, 100, 1.0, 0.03, sig, corr, is_call=False)
    M1 = sum(w[i] * spots[i] * math.exp(0.03 * 1.0) for i in range(3))
    assert abs((c - p) - math.exp(-0.03 * 1.0) * (M1 - 100)) < 1e-9


def test_higher_vol_raises_price():
    spots, w = [100, 90, 110], [0.4, 0.3, 0.3]
    corr = [[1, 0.5, 0.3], [0.5, 1, 0.4], [0.3, 0.4, 1]]
    base = lb(spots, w, 100, 1.0, 0.03, [0.2, 0.25, 0.3], corr)
    high = lb(spots, w, 100, 1.0, 0.03, [0.4, 0.5, 0.6], corr)
    assert high > base


def test_validation():
    with pytest.raises(ValueError):
        lb([100], [1, 2], 105, 1.0, 0.05, [0.2], [[1.0]])   # length mismatch
    with pytest.raises(ValueError):
        lb([100, 90], [0.5, 0.5], 100, 1.0, 0.05, [0.2, 0.2], [[1.0]])  # bad corr shape
