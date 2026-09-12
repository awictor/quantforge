"""Volatility targeting overlay."""

import random

import pytest

from quantforge import target_leverage, vol_targeted_returns, realized_annualized_vol


def test_leverage_formula():
    assert target_leverage(0.10, 0.20) == pytest.approx(0.5)


def test_leverage_capped():
    assert target_leverage(0.10, 0.02, max_leverage=3.0) == 3.0


def test_leverage_monotone_and_floored():
    assert target_leverage(0.10, 0.30) < target_leverage(0.10, 0.15)
    assert target_leverage(0.0, 0.20) == 0.0


def test_overlay_realizes_near_target():
    random.seed(3)
    rets = [random.gauss(0.0005, 0.02) for _ in range(500)]
    overlaid = vol_targeted_returns(rets, 0.15, 60)
    assert abs(realized_annualized_vol(overlaid) - 0.15) < 0.06
    assert realized_annualized_vol(rets) > 0.25   # raw is much higher


def test_validation():
    with pytest.raises(ValueError):
        target_leverage(0.1, 0)
    with pytest.raises(ValueError):
        vol_targeted_returns([0.01, 0.02, 0.03], 0.15, 1)
