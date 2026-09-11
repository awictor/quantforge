"""Cumulative / annualized return and annualized volatility (perfmetrics)."""

import math
import random
import statistics

import pytest

from quantforge import (
    cumulative_return, annualized_return, annualized_volatility,
)


def test_cumulative_return_compounds():
    assert cumulative_return([0.1, -0.1]) == pytest.approx(-0.01, abs=1e-12)


def test_cumulative_return_empty_is_zero():
    assert cumulative_return([]) == 0.0


def test_annualized_return_flat_is_zero():
    assert annualized_return([0.0] * 252) == pytest.approx(0.0, abs=1e-12)


def test_annualized_return_compound():
    r = [0.0004] * 252
    assert annualized_return(r) == pytest.approx(1.0004 ** 252 - 1, abs=1e-9)


def test_annualized_volatility_matches_scaled_stdev():
    rng = random.Random(3)
    r = [rng.gauss(0.0, 0.01) for _ in range(1000)]
    assert annualized_volatility(r) == pytest.approx(
        statistics.stdev(r) * math.sqrt(252), abs=1e-9)


def test_annualized_volatility_flat_is_zero():
    assert annualized_volatility([0.01, 0.01, 0.01]) == pytest.approx(0.0, abs=1e-12)


def test_validation():
    with pytest.raises(ValueError):
        annualized_return([])
    with pytest.raises(ValueError):
        annualized_volatility([0.01])
