"""Continuous-time Markov: generator to transition matrix."""

import pytest

from quantforge import generator_to_transition, generator_default_probability

# 3-state rating generator with an absorbing default state.
Q = [[-0.15, 0.12, 0.03], [0.05, -0.20, 0.15], [0.0, 0.0, 0.0]]


def test_p0_is_identity():
    P = generator_to_transition(Q, 0.0)
    assert all(abs(P[i][j] - (1.0 if i == j else 0.0)) < 1e-12
               for i in range(3) for j in range(3))


def test_rows_sum_to_one_nonneg():
    P = generator_to_transition(Q, 1.0)
    for row in P:
        assert abs(sum(row) - 1.0) < 1e-10
    assert all(v >= 0.0 for row in P for v in row)


def test_semigroup_property():
    P1 = generator_to_transition(Q, 1.0)
    P2 = generator_to_transition(Q, 2.0)
    PP = [[sum(P1[i][k] * P1[k][j] for k in range(3)) for j in range(3)]
          for i in range(3)]
    assert max(abs(PP[i][j] - P2[i][j]) for i in range(3) for j in range(3)) < 1e-9


def test_default_probabilities_increasing():
    dp = generator_default_probability(Q, 2, [1, 2, 5, 10], start_state=0)
    assert all(dp[i] <= dp[i + 1] for i in range(len(dp) - 1))
    assert all(0.0 <= p <= 1.0 for p in dp)


def test_absorbing_state_stays():
    # From the default state the probability of staying is 1 at any horizon.
    P = generator_to_transition(Q, 3.0)
    assert abs(P[2][2] - 1.0) < 1e-10


def test_validation():
    with pytest.raises(ValueError):
        generator_to_transition([[0.1, 0.1], [0.0, 0.0]], 1.0)   # row not summing to 0
    with pytest.raises(ValueError):
        generator_to_transition([[-0.1, -0.05, 0.15], [0, 0, 0], [0, 0, 0]], 1.0)  # neg off-diag
    with pytest.raises(ValueError):
        generator_to_transition(Q, -1.0)
