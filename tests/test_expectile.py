"""Expectile risk measure."""

import random

import pytest

from quantforge import expectile


def _pnl(seed, n=5000):
    rng = random.Random(seed)
    return [rng.gauss(0.001, 0.02) for _ in range(n)]


def test_half_expectile_is_mean_loss():
    pnl = _pnl(1)
    mean_loss = sum(-x for x in pnl) / len(pnl)
    assert abs(expectile(pnl, 0.5) - mean_loss) < 1e-9


def test_monotone_in_tau():
    pnl = _pnl(1)
    es = [expectile(pnl, t) for t in (0.5, 0.9, 0.95, 0.99)]
    assert all(es[i] < es[i + 1] for i in range(len(es) - 1))


def test_first_order_condition():
    pnl = _pnl(2)
    losses = [-x for x in pnl]
    tau = 0.95
    e = expectile(pnl, tau)
    pos = sum(x - e for x in losses if x > e)
    neg = sum(e - x for x in losses if x < e)
    assert abs(tau * pos - (1 - tau) * neg) < 1e-6


def test_above_mean_for_high_tau():
    pnl = _pnl(3)
    mean_loss = sum(-x for x in pnl) / len(pnl)
    assert expectile(pnl, 0.95) > mean_loss


def test_constant_series():
    assert abs(expectile([0.01, 0.01, 0.01], 0.9) - (-0.01)) < 1e-12


def test_validation():
    with pytest.raises(ValueError):
        expectile([], 0.9)
    with pytest.raises(ValueError):
        expectile([0.01, -0.02], 1.0)
