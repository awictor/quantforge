"""Hidden Markov model: forward likelihood and Viterbi decoding."""

import itertools
import math

import pytest

from quantforge import hmm_forward, hmm_viterbi

PI = [0.6, 0.4]
A = [[0.7, 0.3], [0.4, 0.6]]
B = [[0.5, 0.5], [0.1, 0.9]]
OBS = [0, 1, 1, 0]


def _brute_ll(pi, A, B, obs):
    n, T = len(pi), len(obs)
    tot = 0.0
    for path in itertools.product(range(n), repeat=T):
        p = pi[path[0]] * B[path[0]][obs[0]]
        for t in range(1, T):
            p *= A[path[t - 1]][path[t]] * B[path[t]][obs[t]]
        tot += p
    return math.log(tot)


def _brute_viterbi(pi, A, B, obs):
    n, T = len(pi), len(obs)
    best, bp = None, None
    for path in itertools.product(range(n), repeat=T):
        p = pi[path[0]] * B[path[0]][obs[0]]
        for t in range(1, T):
            p *= A[path[t - 1]][path[t]] * B[path[t]][obs[t]]
        if best is None or p > best:
            best, bp = p, path
    return list(bp), math.log(best)


def test_forward_matches_brute_force():
    assert abs(hmm_forward(PI, A, B, OBS) - _brute_ll(PI, A, B, OBS)) < 1e-10


def test_viterbi_matches_brute_force():
    path, lp = hmm_viterbi(PI, A, B, OBS)
    bpath, blp = _brute_viterbi(PI, A, B, OBS)
    assert path == bpath
    assert abs(lp - blp) < 1e-10


def test_deterministic_hmm_recovers_path():
    pi = [1.0, 0.0]
    Ad = [[0.0, 1.0], [1.0, 0.0]]
    Bd = [[1.0, 0.0], [0.0, 1.0]]      # state i emits symbol i
    obs = [0, 1, 0, 1, 0]
    path, _ = hmm_viterbi(pi, Ad, Bd, obs)
    assert path == [0, 1, 0, 1, 0]


def test_forward_scaling_on_long_sequence():
    # A long sequence would underflow an unscaled forward pass; the log-lik stays finite.
    obs = [0, 1] * 500
    ll = hmm_forward(PI, A, B, obs)
    assert math.isfinite(ll)
    assert ll < 0.0


def test_validation():
    with pytest.raises(ValueError):
        hmm_forward(PI, A, B, [])
    with pytest.raises(ValueError):
        hmm_forward(PI, A, B, [0, 5])          # symbol out of range
    with pytest.raises(ValueError):
        hmm_viterbi(PI, [[0.7, 0.3]], B, OBS)  # A not n x n
