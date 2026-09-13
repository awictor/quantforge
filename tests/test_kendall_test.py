"""Tie-corrected Kendall tau-b, Goodman-Kruskal gamma, significance test."""

import math
import random

import pytest

from quantforge import kendall_tau, kendall_tau_b, goodman_kruskal_gamma, kendall_tau_test


def test_perfect_monotone_reaches_one():
    assert abs(kendall_tau_b([1, 2, 3, 4], [5, 6, 7, 8]) - 1.0) < 1e-12
    assert abs(kendall_tau_b([1, 2, 2, 3], [10, 20, 20, 30]) - 1.0) < 1e-12   # aligned ties
    assert abs(kendall_tau_b([1, 2, 3, 4], [8, 7, 6, 5]) + 1.0) < 1e-12


def test_tau_b_equals_tau_a_without_ties():
    rng = random.Random(1)
    for _ in range(300):
        n = rng.randint(2, 15)
        xs = random.Random(rng.random()).sample(range(1000), n)
        ys = random.Random(rng.random()).sample(range(1000), n)
        assert abs(kendall_tau(xs, ys) - kendall_tau_b(xs, ys)) < 1e-9


def test_gamma_formula():
    # x=[1,2,3,4] y=[1,2,4,3]: C=5, D=1 -> gamma = 4/6.
    assert abs(goodman_kruskal_gamma([1, 2, 3, 4], [1, 2, 4, 3]) - 2.0 / 3.0) < 1e-12


def test_all_tied_margin_is_zero():
    assert kendall_tau_b([5, 5, 5], [1, 2, 3]) == 0.0


def test_significance_independence_vs_dependence():
    rng = random.Random(5)
    xs = [rng.gauss(0, 1) for _ in range(60)]
    ys = [rng.gauss(0, 1) for _ in range(60)]
    r = kendall_tau_test(xs, ys)
    assert r["p_value"] > 0.05          # independent -> not significant

    mono = list(range(60))
    r2 = kendall_tau_test(mono, [2 * v for v in mono])
    assert abs(r2["tau_b"] - 1.0) < 1e-12
    assert r2["p_value"] < 1e-10        # perfect dependence -> highly significant


def test_p_value_matches_erfc():
    rng = random.Random(9)
    xs = [rng.gauss(0, 1) for _ in range(40)]
    ys = [rng.gauss(0, 1) for _ in range(40)]
    r = kendall_tau_test(xs, ys)
    assert abs(r["p_value"] - math.erfc(abs(r["z"]) / math.sqrt(2.0))) < 1e-12


def test_variance_formula():
    # Var(S) = n(n-1)(2n+5)/18; for a perfectly monotone n=10, S = 45.
    r = kendall_tau_test(list(range(10)), list(range(10)))
    var = 10 * 9 * 25 / 18.0
    assert abs(r["z"] - 45.0 / math.sqrt(var)) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        kendall_tau_b([1], [1])
    with pytest.raises(ValueError):
        kendall_tau_b([1, 2], [1])
