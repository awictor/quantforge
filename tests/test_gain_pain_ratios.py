"""Gain-to-pain, Sterling, and Burke ratios."""

import math

import pytest

from quantforge import gain_to_pain_ratio, sterling_ratio, burke_ratio
from quantforge.perfmetrics import drawdown_curve

R = [0.02, -0.01, 0.03, -0.02, 0.01, -0.015, 0.025, 0.005, -0.008, 0.012]


def test_gain_to_pain_matches_manual():
    total = sum(R)
    loss = sum(-x for x in R if x < 0)
    assert abs(gain_to_pain_ratio(R) - total / loss) < 1e-12


def test_gain_to_pain_infinite_without_losses():
    assert gain_to_pain_ratio([0.01, 0.02, 0.03]) == float("inf")


def test_sterling_matches_manual():
    dd = drawdown_curve(R)
    avg_dd = sum(dd) / len(dd)
    ann = (sum(R) / len(R)) * 252
    assert abs(sterling_ratio(R) - ann / (avg_dd + 0.10)) < 1e-9


def test_burke_matches_manual():
    dd = drawdown_curve(R)
    ss = math.sqrt(sum(d * d for d in dd))
    ann = (sum(R) / len(R)) * 252
    assert abs(burke_ratio(R) - ann / ss) < 1e-9


def test_burke_penalizes_deep_drawdowns():
    good = [0.01] * 10 + [-0.05] + [0.01] * 10
    bad = [0.01] * 10 + [-0.30] + [0.01] * 10
    assert burke_ratio(good) > burke_ratio(bad)


def test_validation():
    with pytest.raises(ValueError):
        gain_to_pain_ratio([])
    with pytest.raises(ValueError):
        burke_ratio([0.01, 0.02, 0.03])       # no drawdown
