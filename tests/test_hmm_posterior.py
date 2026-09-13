"""HMM forward-backward posterior state probabilities."""

import itertools

import pytest

from quantforge import hmm_posterior

PI = [0.6, 0.4]
A = [[0.7, 0.3], [0.4, 0.6]]
B = [[0.5, 0.5], [0.1, 0.9]]
OBS = [0, 1, 1, 0]


def _brute_gamma(pi, A, B, obs, t, i):
    n, T = len(pi), len(obs)
    num = den = 0.0
    for path in itertools.product(range(n), repeat=T):
        p = pi[path[0]] * B[path[0]][obs[0]]
        for s in range(1, T):
            p *= A[path[s - 1]][path[s]] * B[path[s]][obs[s]]
        den += p
        if path[t] == i:
            num += p
    return num / den


def test_rows_sum_to_one():
    g = hmm_posterior(PI, A, B, OBS)
    assert all(abs(sum(row) - 1.0) < 1e-10 for row in g)


def test_matches_brute_force_marginals():
    g = hmm_posterior(PI, A, B, OBS)
    for t in range(len(OBS)):
        for i in range(2):
            assert abs(g[t][i] - _brute_gamma(PI, A, B, OBS, t, i)) < 1e-9


def test_deterministic_concentrates():
    pi = [1.0, 0.0]
    Ad = [[0.0, 1.0], [1.0, 0.0]]
    Bd = [[1.0, 0.0], [0.0, 1.0]]
    g = hmm_posterior(pi, Ad, Bd, [0, 1, 0, 1])
    assert g == [[1.0, 0.0], [0.0, 1.0], [1.0, 0.0], [0.0, 1.0]]


def test_length_matches_observations():
    g = hmm_posterior(PI, A, B, [0, 1, 0])
    assert len(g) == 3


def test_validation():
    with pytest.raises(ValueError):
        hmm_posterior(PI, A, B, [])
    with pytest.raises(ValueError):
        hmm_posterior(PI, A, B, [0, 9])
