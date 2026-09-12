"""Jump-robust realized volatility measures."""

import math
import random

import pytest

from quantforge import (
    realized_variance_from_returns, bipower_variation, jump_variation,
    realized_volatility_signature,
)


def _continuous(n, sigma, seed):
    rng = random.Random(seed)
    return [rng.gauss(0, sigma) for _ in range(n)]


def test_rv_and_bv_estimate_integrated_variance():
    sigma, n = 0.02, 2000
    ret = _continuous(n, sigma, 1)
    iv = n * sigma * sigma
    assert abs(realized_variance_from_returns(ret) - iv) / iv < 0.1
    assert abs(bipower_variation(ret) - iv) / iv < 0.1


def test_continuous_path_has_negligible_jump():
    ret = _continuous(2000, 0.02, 1)
    rv = realized_variance_from_returns(ret)
    assert jump_variation(ret) / rv < 0.1


def test_jump_inflates_rv_but_not_bv():
    ret = _continuous(2000, 0.02, 1)
    retj = list(ret)
    retj[1000] += 0.5                      # inject a jump
    rv0 = realized_variance_from_returns(ret)
    rvj = realized_variance_from_returns(retj)
    bv0 = bipower_variation(ret)
    bvj = bipower_variation(retj)
    assert rvj > rv0 + 0.2                  # RV absorbs the jump^2 = 0.25
    assert abs(bvj - bv0) / bv0 < 0.15      # BV is jump-robust


def test_jump_variation_positive_and_near_jump_squared():
    ret = _continuous(2000, 0.02, 1)
    retj = list(ret)
    retj[1000] += 0.5
    jv = jump_variation(retj)
    assert jv > 0.15
    assert abs(jv - 0.25) / 0.25 < 0.25     # ~ jump^2, minus a little BV leakage


def test_signature_is_sqrt_of_rv():
    ret = _continuous(500, 0.02, 3)
    rv = realized_variance_from_returns(ret)
    assert abs(realized_volatility_signature(ret) - math.sqrt(rv)) < 1e-12
    assert abs(realized_volatility_signature(ret, 252) - math.sqrt(rv * 252)) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        realized_variance_from_returns([])
    with pytest.raises(ValueError):
        bipower_variation([0.1])            # need >= 2
