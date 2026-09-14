"""Tests for TwoSat, cross-checked against brute-force assignment enumeration."""

import itertools
import random

import pytest

from quantforge.two_sat import TwoSat


def _brute_sat(n, clauses):
    for bits in itertools.product([False, True], repeat=n):
        if all(
            (bits[a] if a >= 0 else not bits[~a]) or (bits[b] if b >= 0 else not bits[~b])
            for a, b in clauses
        ):
            return list(bits)
    return None


def _check(clauses, assign):
    return all(
        (assign[a] if a >= 0 else not assign[~a]) or (assign[b] if b >= 0 else not assign[~b])
        for a, b in clauses
    )


def test_fuzz_satisfiability_agrees_and_assignment_valid():
    rng = random.Random(71)
    sat = unsat = 0
    for _ in range(5000):
        n = rng.randint(1, 6)
        m = rng.randint(0, 12)
        clauses = []
        for _ in range(m):
            a = rng.randint(0, n - 1)
            b = rng.randint(0, n - 1)
            if rng.random() < 0.5:
                a = ~a
            if rng.random() < 0.5:
                b = ~b
            clauses.append((a, b))
        ts = TwoSat(n)
        for a, b in clauses:
            ts.add_or(a, b)
        res = ts.solve()
        bru = _brute_sat(n, clauses)
        assert (res is None) == (bru is None)
        if res is not None:
            assert _check(clauses, res)
            sat += 1
        else:
            unsat += 1
    assert sat > 0 and unsat > 0  # the fuzz exercised both outcomes


def test_contradiction_is_unsatisfiable():
    ts = TwoSat(1)
    ts.force_true(0)
    ts.force_true(~0)
    assert ts.solve() is None
    assert ts.is_satisfiable() is False


def test_simple_satisfiable():
    clauses = [(0, 1), (~0, 2), (~1, ~2)]
    ts = TwoSat(3)
    for a, b in clauses:
        ts.add_or(a, b)
    res = ts.solve()
    assert res is not None
    assert _check(clauses, res)


def test_force_true_sets_variable():
    ts = TwoSat(2)
    ts.force_true(0)
    ts.force_true(~1)  # force var 1 false
    res = ts.solve()
    assert res[0] is True
    assert res[1] is False


def test_implication_chain():
    ts = TwoSat(3)
    ts.force_true(0)
    ts.add_implication(0, 1)  # 0 -> 1
    ts.add_implication(1, 2)  # 1 -> 2
    res = ts.solve()
    assert res == [True, True, True]


def test_no_clauses_any_assignment_valid():
    res = TwoSat(3).solve()
    assert isinstance(res, list)
    assert len(res) == 3
    assert all(isinstance(b, bool) for b in res)


def test_zero_variables():
    assert TwoSat(0).solve() == []
    assert TwoSat(0).is_satisfiable() is True


def test_is_satisfiable_matches_solve():
    ts = TwoSat(2)
    ts.add_or(0, 1)
    assert ts.is_satisfiable() is (ts.solve() is not None)


def test_negative_n_raises():
    with pytest.raises(ValueError):
        TwoSat(-1)


def test_out_of_range_literal_raises():
    ts = TwoSat(2)
    with pytest.raises(IndexError):
        ts.add_or(0, 5)
    with pytest.raises(IndexError):
        ts.add_or(~9, 0)


def test_add_or_returns_self_for_chaining():
    ts = TwoSat(2)
    assert ts.add_or(0, 1) is ts


def test_both_negated_clause():
    # (~0 OR ~1): forbids 0=1=True
    ts = TwoSat(2)
    ts.force_true(0)
    ts.force_true(1)
    ts.add_or(~0, ~1)
    assert ts.solve() is None
