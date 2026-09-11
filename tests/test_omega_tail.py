"""Omega ratio and tail ratio (perfmetrics module)."""

import pytest

from quantforge import omega_ratio, tail_ratio


R = [0.01, -0.005, 0.02, 0.008, -0.012, 0.015, 0.003, -0.002]


def test_omega_matches_manual():
    up = sum(max(x, 0.0) for x in R)
    down = sum(max(-x, 0.0) for x in R)
    assert omega_ratio(R) == pytest.approx(up / down, abs=1e-12)


def test_omega_above_one_for_net_positive():
    assert omega_ratio(R) > 1.0


def test_omega_infinite_without_downside():
    assert omega_ratio([0.01, 0.02]) == float("inf")


def test_omega_symmetric_is_one():
    assert omega_ratio([0.01, -0.01, 0.02, -0.02]) == pytest.approx(1.0, abs=1e-12)


def test_tail_ratio_symmetric_is_one():
    assert tail_ratio([0.01, -0.01, 0.02, -0.02]) == pytest.approx(1.0, abs=1e-9)


def test_tail_ratio_right_skew_above_one():
    skew = [0.05, -0.01, -0.01, -0.01, 0.01, 0.02]
    assert tail_ratio(skew, pct=20.0) > 1.0


def test_validation():
    with pytest.raises(ValueError):
        omega_ratio([0.0, 0.0])       # all at threshold
    with pytest.raises(ValueError):
        tail_ratio(R, pct=60.0)       # pct out of range
    with pytest.raises(ValueError):
        omega_ratio([])
