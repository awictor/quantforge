"""HMM Baum-Welch parameter estimation and simulation."""

import pytest

from quantforge import hmm_simulate, hmm_baum_welch, hmm_forward

PI = [0.5, 0.5]
A = [[0.9, 0.1], [0.1, 0.9]]
B = [[0.9, 0.1], [0.1, 0.9]]


def _obs():
    _, obs = hmm_simulate(PI, A, B, 4000, seed=42)
    return obs


def test_recovers_parameters_up_to_permutation():
    obs = _obs()
    r = hmm_baum_welch(obs, 2, 2, max_iter=100, seed=7)
    Be = r["B"]
    s0 = 0 if Be[0][0] > Be[1][0] else 1
    s1 = 1 - s0
    assert abs(Be[s0][0] - 0.9) < 0.06
    assert abs(Be[s1][1] - 0.9) < 0.06
    assert abs(r["A"][s0][s0] - 0.9) < 0.06
    assert abs(r["A"][s1][s1] - 0.9) < 0.06


def test_fit_likelihood_beats_true_model():
    obs = _obs()
    r = hmm_baum_welch(obs, 2, 2, max_iter=100, seed=7)
    assert r["log_likelihood"] >= hmm_forward(PI, A, B, obs) - 2.0


def test_log_likelihood_non_decreasing():
    obs = _obs()
    lls = [hmm_baum_welch(obs, 2, 2, max_iter=m, seed=7)["log_likelihood"]
           for m in (1, 3, 5, 10, 30)]
    assert all(lls[i] <= lls[i + 1] + 1e-4 for i in range(len(lls) - 1))


def test_simulate_shapes():
    states, obs = hmm_simulate(PI, A, B, 100, seed=1)
    assert len(states) == len(obs) == 100
    assert all(s in (0, 1) for s in states)
    assert all(o in (0, 1) for o in obs)


def test_validation():
    with pytest.raises(ValueError):
        hmm_baum_welch([0], 2, 2)
    with pytest.raises(ValueError):
        hmm_baum_welch([0, 5], 2, 2)      # symbol out of range
