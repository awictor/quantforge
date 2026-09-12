"""Finite-state Markov chains."""

import pytest

from quantforge import (
    n_step_transition, stationary_distribution, expected_hitting_time,
)


P = [[0.9, 0.1], [0.2, 0.8]]


def test_stationary_two_state_closed_form():
    pi = stationary_distribution(P)
    assert pi[0] == pytest.approx(2 / 3, abs=1e-6)
    assert pi[1] == pytest.approx(1 / 3, abs=1e-6)


def test_stationary_fixed_point():
    pi = stationary_distribution(P)
    piP = [sum(pi[i] * P[i][j] for i in range(2)) for j in range(2)]
    assert all(piP[j] == pytest.approx(pi[j], abs=1e-9) for j in range(2))
    assert sum(pi) == pytest.approx(1.0)


def test_n_step_row_stochastic():
    P5 = n_step_transition(P, 5)
    assert all(sum(row) == pytest.approx(1.0) for row in P5)


def test_n_step_identity_and_one():
    assert n_step_transition(P, 0) == [[1.0, 0.0], [0.0, 1.0]]
    assert n_step_transition(P, 1) == P


def test_n_step_converges_to_stationary():
    pi = stationary_distribution(P)
    P100 = n_step_transition(P, 100)
    assert all(P100[0][j] == pytest.approx(pi[j], abs=1e-6) for j in range(2))


def test_expected_hitting_time():
    h = expected_hitting_time(P, 0)
    assert h[0] == 0.0
    assert h[1] == pytest.approx(5.0, abs=1e-6)  # 1 + 0.8 h1 => h1 = 5


def test_validation():
    with pytest.raises(ValueError):
        stationary_distribution([[0.5, 0.6], [0.2, 0.8]])   # row not summing to 1
    with pytest.raises(ValueError):
        n_step_transition(P, -1)
    with pytest.raises(ValueError):
        expected_hitting_time(P, 5)
