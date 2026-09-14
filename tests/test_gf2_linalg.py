"""Tests for GF(2) linear algebra, cross-checked against brute-force enumeration."""

import itertools
import math
import random

import pytest

from quantforge.gf2_linalg import solve_gf2, gf2_rank, gf2_nullspace_basis


def _eval_eq(mask, x):
    v = 0
    j = 0
    m = mask
    while m:
        if m & 1:
            v ^= x[j]
        m >>= 1
        j += 1
    return v


def _brute_has_solution(eqs, rhs, n):
    for bits in itertools.product([0, 1], repeat=n):
        if all(_eval_eq(eqs[i], bits) == rhs[i] for i in range(len(eqs))):
            return True
    return False


def _brute_rank(rows):
    span = {0}
    for r in rows:
        span |= {s ^ r for s in span}
    return int(round(math.log2(len(span)))) if len(span) > 1 else 0


def test_fuzz_solve_and_rank_and_nullspace():
    rng = random.Random(211)
    for _ in range(4000):
        n = rng.randint(1, 7)
        m = rng.randint(0, 8)
        eqs = [rng.randint(0, (1 << n) - 1) for _ in range(m)]
        rhs = [rng.randint(0, 1) for _ in range(m)]

        sol = solve_gf2(eqs, rhs, n)
        assert (sol is not None) == _brute_has_solution(eqs, rhs, n)
        if sol is not None:
            assert all(_eval_eq(eqs[i], sol) == rhs[i] for i in range(m))

        rank = gf2_rank(eqs)
        assert rank == _brute_rank(eqs)

        nb = gf2_nullspace_basis(eqs, n)
        assert len(nb) == n - rank
        for vec in nb:
            assert all(_eval_eq(eqs[i], vec) == 0 for i in range(m))
        if nb:
            vints = [sum(v[j] << j for j in range(n)) for v in nb]
            assert gf2_rank(vints) == len(nb)  # basis vectors independent


def test_explicit_unique_solution():
    eqs = [0b011, 0b110, 0b001]  # x0^x1=1, x1^x2=0, x0=1
    rhs = [1, 0, 1]
    assert solve_gf2(eqs, rhs, 3) == [1, 0, 0]


def test_inconsistent_returns_none():
    assert solve_gf2([0b1, 0b1], [0, 1], 1) is None


def test_underdetermined_returns_a_valid_solution():
    s = solve_gf2([0b11], [1], 2)  # x0 ^ x1 = 1
    assert _eval_eq(0b11, s) == 1


def test_rank_full():
    assert gf2_rank([0b001, 0b010, 0b100]) == 3
    assert gf2_rank([0b11, 0b11]) == 1  # duplicate rows


def test_rank_empty():
    assert gf2_rank([]) == 0
    assert gf2_rank([0, 0]) == 0


def test_nullspace_empty_system_is_full_space():
    nb = gf2_nullspace_basis([], 3)
    assert len(nb) == 3


def test_nullspace_full_rank_is_trivial():
    nb = gf2_nullspace_basis([0b001, 0b010, 0b100], 3)
    assert nb == []


def test_nullspace_vectors_satisfy_system():
    eqs = [0b0011, 0b1100]  # x0^x1=0, x2^x3=0
    nb = gf2_nullspace_basis(eqs, 4)
    assert len(nb) == 2
    for vec in nb:
        assert all(_eval_eq(e, vec) == 0 for e in eqs)


def test_length_mismatch_raises():
    with pytest.raises(ValueError):
        solve_gf2([1, 2], [1], 2)


def test_zero_solution_when_rhs_all_zero():
    sol = solve_gf2([0b011, 0b110], [0, 0], 3)
    assert sol is not None
    assert all(_eval_eq(e, sol) == 0 for e in (0b011, 0b110))


def test_single_variable():
    assert solve_gf2([0b1], [1], 1) == [1]
    assert solve_gf2([0b1], [0], 1) == [0]
