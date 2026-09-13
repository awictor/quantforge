"""Portfolio concentration and diversification measures."""

import random

import pytest

from quantforge import (herfindahl_index, effective_number_of_constituents,
                        effective_number_of_bets)


def test_equal_weight_hhi_and_enc():
    w = [0.25] * 4
    assert abs(herfindahl_index(w) - 0.25) < 1e-12
    assert abs(effective_number_of_constituents(w) - 4.0) < 1e-9


def test_concentration_approaches_one():
    enc = effective_number_of_constituents([0.97, 0.01, 0.01, 0.01])
    assert enc < 1.2


def test_enb_identity_equals_n():
    w = [0.25] * 4
    I = [[1.0 if i == j else 0.0 for j in range(4)] for i in range(4)]
    assert abs(effective_number_of_bets(w, I) - 4.0) < 1e-6


def test_enb_dominant_factor_near_one():
    cov = [[100.0 if i == j == 0 else (1.0 if i == j else 0.0) for j in range(4)]
           for i in range(4)]
    assert effective_number_of_bets([0.25] * 4, cov) < 1.5


def test_enb_bounded_by_n():
    rng = random.Random(1)
    A = [[rng.gauss(0, 1) for _ in range(4)] for _ in range(4)]
    cov = [[sum(A[i][k] * A[j][k] for k in range(4)) for j in range(4)]
           for i in range(4)]
    enb = effective_number_of_bets([0.25] * 4, cov)
    assert 1.0 <= enb <= 4.0 + 1e-9


def test_hhi_normalizes_weights():
    # Unnormalized weights give the same HHI as their normalized form.
    assert abs(herfindahl_index([2, 2, 2, 2]) - 0.25) < 1e-12


def test_validation():
    with pytest.raises(ValueError):
        herfindahl_index([])
    with pytest.raises(ValueError):
        herfindahl_index([0.0, 0.0])
    with pytest.raises(ValueError):
        effective_number_of_bets([0.5, 0.5], [[1.0, 0.0]])   # wrong cov shape
