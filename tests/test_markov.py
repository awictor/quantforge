"""Finite-state Markov chains."""

import pytest

from quantforge import (
    n_step_transition, stationary_distribution, expected_hitting_time,
    fundamental_matrix, expected_steps_to_absorption, absorption_probabilities,
    cumulative_default_term_structure, marginal_default_probabilities,
)


RATING_P = [[0.90, 0.09, 0.01], [0.05, 0.85, 0.10], [0.0, 0.0, 1.0]]
HZ = [1, 2, 3, 5, 10]


ABS_P = [[0.5, 0.3, 0.2], [0.1, 0.6, 0.3], [0.0, 0.0, 1.0]]
TRANS = [0, 1]


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


def test_cumulative_default_monotone_and_matches_nstep():
    cum = cumulative_default_term_structure(RATING_P, 2, HZ)
    assert all(cum[i] <= cum[i + 1] for i in range(len(cum) - 1))
    assert all(cum[i] == pytest.approx(n_step_transition(RATING_P, HZ[i])[0][2])
               for i in range(len(HZ)))
    assert cum[0] == pytest.approx(0.01)


def test_marginal_defaults_nonneg_and_sum():
    cum = cumulative_default_term_structure(RATING_P, 2, HZ)
    marg = marginal_default_probabilities(RATING_P, 2, HZ)
    assert all(m >= 0 for m in marg)
    assert sum(marg) == pytest.approx(cum[-1])


def test_riskier_start_higher_default():
    cumA = cumulative_default_term_structure(RATING_P, 2, HZ, start_state=0)
    cumB = cumulative_default_term_structure(RATING_P, 2, HZ, start_state=1)
    assert all(cumB[i] >= cumA[i] for i in range(len(HZ)))


def test_default_converges_to_one():
    assert cumulative_default_term_structure(RATING_P, 2, [200])[0] > 0.99


def test_default_requires_absorbing():
    with pytest.raises(ValueError):
        cumulative_default_term_structure([[0.9, 0.1], [0.1, 0.9]], 1, [1])


def test_fundamental_matrix_inverts_i_minus_q():
    N = fundamental_matrix(ABS_P, TRANS)
    Q = [[ABS_P[i][j] for j in TRANS] for i in TRANS]
    prod = [[sum(((1 if r == c else 0) - Q[r][c]) * N[c][cc] for c in range(2))
             for cc in range(2)] for r in range(2)]
    assert all(prod[i][j] == pytest.approx(1.0 if i == j else 0.0, abs=1e-9)
               for i in range(2) for j in range(2))


def test_expected_steps_matches_hitting_time():
    steps = expected_steps_to_absorption(ABS_P, TRANS)
    h = expected_hitting_time(ABS_P, 2)
    assert steps[0] == pytest.approx(h[0], abs=1e-6)
    assert all(s > 0 for s in steps)


def test_absorption_probabilities_sum_to_one():
    B = absorption_probabilities(ABS_P, TRANS, [2])
    assert all(sum(row) == pytest.approx(1.0) for row in B)


def test_two_absorbing_states():
    P2 = [[0.4, 0.3, 0.2, 0.1], [0.2, 0.5, 0.1, 0.2], [0, 0, 1, 0], [0, 0, 0, 1]]
    B = absorption_probabilities(P2, [0, 1], [2, 3])
    assert all(sum(row) == pytest.approx(1.0) for row in B)


def test_absorbing_validation():
    with pytest.raises(ValueError):
        fundamental_matrix(ABS_P, [])


def test_validation():
    with pytest.raises(ValueError):
        stationary_distribution([[0.5, 0.6], [0.2, 0.8]])   # row not summing to 1
    with pytest.raises(ValueError):
        n_step_transition(P, -1)
    with pytest.raises(ValueError):
        expected_hitting_time(P, 5)
