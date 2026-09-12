"""Tail-factor extrapolation for chain-ladder."""

import pytest

from quantforge import exponential_tail_factor, chain_ladder_with_tail
from quantforge.chain_ladder import chain_ladder


TRI = [[100, 150, 180], [110, 165], [120]]


def test_decaying_factors_finite_tail_above_one():
    t = exponential_tail_factor([1.5, 1.2, 1.1, 1.05])
    assert 1.0 < t < 2.0


def test_flat_factors_tail_one():
    assert exponential_tail_factor([1.0, 1.0]) == 1.0


def test_tail_raises_ultimates_proportionally():
    base = chain_ladder(TRI)
    wt = chain_ladder_with_tail(TRI, 1.05)
    assert all(abs(wt["ultimate"][i] - base["ultimate"][i] * 1.05) < 1e-9
               for i in range(len(TRI)))


def test_unit_tail_reproduces_chain_ladder():
    base = chain_ladder(TRI)
    wt = chain_ladder_with_tail(TRI, 1.0)
    assert all(abs(wt["ultimate"][i] - base["ultimate"][i]) < 1e-9
               for i in range(len(TRI)))


def test_validation():
    with pytest.raises(ValueError):
        exponential_tail_factor([1.1, 1.2, 1.3])   # non-decaying excess
    with pytest.raises(ValueError):
        chain_ladder_with_tail(TRI, 0.9)            # tail < 1
