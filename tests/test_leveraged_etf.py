"""Leveraged / inverse ETF path and volatility drag."""

import pytest

from quantforge import (
    leveraged_etf_path, volatility_drag, expected_leveraged_return,
    flat_market_decay,
)


def test_one_x_matches_underlying():
    rets = [0.02, -0.01, 0.03, -0.02]
    und = 1.0
    for r in rets:
        und *= (1 + r)
    assert leveraged_etf_path(rets, 1.0)[-1] == pytest.approx(und)


def test_drag_formula_and_zeros():
    assert volatility_drag(2, 0.20) == pytest.approx(0.5 * 2 * 1 * 0.04)
    assert volatility_drag(1, 0.2) == 0
    assert volatility_drag(0, 0.2) == 0


def test_inverse_and_higher_leverage_more_drag():
    assert volatility_drag(-1, 0.2) > 0
    assert volatility_drag(3, 0.2) > volatility_drag(2, 0.2)


def test_expected_return_below_naive():
    mu = 0.0004
    assert expected_leveraged_return(mu, 3, 0.02, 252) < 3 * mu * 252


def test_flat_market_decay():
    rt = [0.10, 1 / 1.1 - 1]   # exact round trip
    assert flat_market_decay(3, rt) < 0
    assert flat_market_decay(1, rt) == pytest.approx(0.0, abs=1e-9)


def test_validation():
    with pytest.raises(ValueError):
        volatility_drag(2, -1)
    with pytest.raises(ValueError):
        expected_leveraged_return(0.001, 2, 0.2, -1)
