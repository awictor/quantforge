"""Stable log-sum-exp, softmax, log-softmax."""

import math

import pytest

from quantforge import logsumexp, softmax, log_softmax


def test_matches_naive_small():
    x = [1.0, 2.0, 3.0]
    naive = math.log(sum(math.exp(v) for v in x))
    assert abs(logsumexp(x) - naive) < 1e-12


def test_no_overflow_large():
    assert abs(logsumexp([1000.0, 1001.0, 1002.0]) - 1002.4076059644444) < 1e-6


def test_no_underflow_tiny():
    assert logsumexp([-1000.0, -1001.0]) < -900          # finite, not -inf


def test_softmax_sums_to_one():
    s = softmax([1.0, 2.0, 3.0])
    assert abs(sum(s) - 1.0) < 1e-12
    e = [math.exp(v) for v in [1, 2, 3]]
    tot = sum(e)
    assert max(abs(s[i] - e[i] / tot) for i in range(3)) < 1e-12


def test_softmax_stable_large():
    s = softmax([1000.0, 1001.0, 1002.0])
    assert abs(sum(s) - 1.0) < 1e-10


def test_log_softmax_consistency():
    ls = log_softmax([1.0, 2.0, 3.0])
    s = softmax([1.0, 2.0, 3.0])
    assert max(abs(ls[i] - math.log(s[i])) for i in range(3)) < 1e-12
    assert abs(sum(math.exp(v) for v in ls) - 1.0) < 1e-12


def test_weighted():
    w = [2.0, 1.0, 3.0]
    x = [0.0, 1.0, 2.0]
    naive = math.log(sum(w[i] * math.exp(x[i]) for i in range(3)))
    assert abs(logsumexp(x, weights=w) - naive) < 1e-12


def test_uniform_softmax():
    assert softmax([5.0, 5.0, 5.0]) == [1 / 3, 1 / 3, 1 / 3]


def test_validation():
    with pytest.raises(ValueError):
        logsumexp([])
    with pytest.raises(ValueError):
        logsumexp([1.0, 2.0], weights=[1.0])
    with pytest.raises(ValueError):
        logsumexp([1.0], weights=[-1.0])
