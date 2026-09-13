"""Drawdown risk measures: DaR, CDaR, average drawdown."""

import random

import pytest

from quantforge import (average_drawdown, drawdown_at_risk,
                        conditional_drawdown_at_risk)
from quantforge.drawdown import drawdown_analytics


def _rets(seed=1, n=2000):
    rng = random.Random(seed)
    return [rng.gauss(0.0005, 0.02) for _ in range(n)]


def test_ordering():
    r = _rets()
    avg = average_drawdown(r)
    dar = drawdown_at_risk(r, 0.95)
    cdar = conditional_drawdown_at_risk(r, 0.95)
    assert cdar >= dar >= avg >= 0.0


def test_bounded_by_max_drawdown():
    r = _rets()
    maxdd = drawdown_analytics(r)["max_drawdown_depth"]
    assert conditional_drawdown_at_risk(r, 0.95) <= maxdd + 1e-9
    assert drawdown_at_risk(r, 0.95) <= maxdd + 1e-9


def test_monotone_in_confidence():
    r = _rets()
    assert drawdown_at_risk(r, 0.90) <= drawdown_at_risk(r, 0.95) <= drawdown_at_risk(r, 0.99)
    assert (conditional_drawdown_at_risk(r, 0.90)
            <= conditional_drawdown_at_risk(r, 0.99))


def test_monotonic_decline_equals_max():
    dec = [-0.01] * 100
    cdar = conditional_drawdown_at_risk(dec, 0.99)
    maxdd = drawdown_analytics(dec)["max_drawdown_depth"]
    assert abs(cdar - maxdd) < 1e-9


def test_no_drawdown_when_monotone_up():
    up = [0.01] * 50
    assert average_drawdown(up) == 0.0
    assert conditional_drawdown_at_risk(up, 0.95) == 0.0


def test_validation():
    with pytest.raises(ValueError):
        average_drawdown([])
    with pytest.raises(ValueError):
        drawdown_at_risk([0.01, -0.02], confidence=1.0)
    with pytest.raises(ValueError):
        conditional_drawdown_at_risk([0.01, -0.02], confidence=0.0)
