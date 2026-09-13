"""Detrended fluctuation analysis."""

import random

import pytest

from quantforge import dfa_exponent, dfa_fluctuations


def _white(n, seed):
    rng = random.Random(seed)
    return [rng.gauss(0, 1) for _ in range(n)]


def _walk(n, seed):
    rng = random.Random(seed)
    out, acc = [], 0.0
    for _ in range(n):
        acc += rng.gauss(0, 1)
        out.append(acc)
    return out


def test_white_noise_half():
    a = [dfa_exponent(_white(4000, s)) for s in range(12)]
    assert abs(sum(a) / len(a) - 0.5) < 0.06


def test_random_walk_one_and_a_half():
    a = [dfa_exponent(_walk(4000, s)) for s in range(12)]
    assert abs(sum(a) / len(a) - 1.5) < 0.1


def test_anti_persistent_below_half():
    a = []
    for s in range(12):
        wn = _white(4000, s)
        diff = [wn[i] - wn[i - 1] for i in range(1, len(wn))]
        a.append(dfa_exponent(diff))
    assert sum(a) / len(a) < 0.4


def test_integration_adds_one():
    # DFA(cumsum(x)) ~ DFA(x) + 1.
    wn = _white(4000, 3)
    walk = []
    acc = 0.0
    for v in wn:
        acc += v
        walk.append(acc)
    assert dfa_exponent(walk) - dfa_exponent(wn) > 0.8


def test_fluctuations_increase_with_scale():
    scales, flucts = dfa_fluctuations(_walk(2000, 1))
    # For a persistent series the fluctuation grows with the window size.
    assert flucts[-1] > flucts[0]
    assert len(scales) == len(flucts) >= 2


def test_validation():
    with pytest.raises(ValueError):
        dfa_exponent([1.0, 2.0, 3.0])          # < 16 points
    with pytest.raises(ValueError):
        dfa_fluctuations(list(range(10)))
