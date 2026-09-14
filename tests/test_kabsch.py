"""Tests for the Kabsch alignment, cross-checked by recovering known rigid transforms."""

import math
import random

import pytest

from quantforge.kabsch import kabsch
from quantforge.so3 import rodrigues


def _close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def _apply(R, t, p):
    return [sum(R[k][j] * p[j] for j in range(3)) + t[k] for k in range(3)]


def _det3(R):
    return (R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1])
            - R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0])
            + R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]))


def test_fuzz_exact_recovery():
    rng = random.Random(471)
    for _ in range(3000):
        n = rng.randint(3, 15)
        P = [[rng.uniform(-10, 10) for _ in range(3)] for _ in range(n)]
        omega = [rng.uniform(-math.pi, math.pi) for _ in range(3)]
        Rt = rodrigues(omega)
        tt = [rng.uniform(-5, 5) for _ in range(3)]
        Q = [_apply(Rt, tt, p) for p in P]
        R, t, rmsd = kabsch(P, Q)
        assert rmsd < 1e-6
        for i in range(n):
            pred = _apply(R, t, P[i])
            assert all(_close(pred[k], Q[i][k], 1e-5) for k in range(3))
        for a in range(3):
            for b in range(3):
                d = sum(R[a][k] * R[b][k] for k in range(3))
                assert _close(d, 1.0 if a == b else 0.0, 1e-6)
        assert _close(_det3(R), 1.0, 1e-6)


def test_noisy_alignment_small_rmsd():
    rng = random.Random(472)
    P = [[rng.uniform(-5, 5) for _ in range(3)] for _ in range(20)]
    Rt = rodrigues([0.3, 0.5, -0.2])
    tt = [1, 2, 3]
    Q = [[_apply(Rt, tt, p)[k] + rng.gauss(0, 0.01) for k in range(3)] for p in P]
    R, t, rmsd = kabsch(P, Q)
    assert 0 < rmsd < 0.1


def test_identity_when_equal():
    P = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]
    R, t, rmsd = kabsch(P, P)
    assert rmsd < 1e-9
    for i in range(3):
        assert _close(R[i][i], 1.0)


def test_pure_translation():
    P = [[1, 0, 0], [0, 2, 0], [0, 0, 3], [1, 1, 1], [2, -1, 0]]  # non-degenerate
    shift = [10, -5, 2]
    Q = [[p[k] + shift[k] for k in range(3)] for p in P]
    R, t, rmsd = kabsch(P, Q)
    assert rmsd < 1e-9
    # non-degenerate + no rotation -> R is identity and t is the shift
    for i in range(3):
        assert _close(R[i][i], 1.0)
    assert all(_close(t[k], shift[k], 1e-6) for k in range(3))
    assert _close(_det3(R), 1.0)


def test_recovers_90_degree_rotation():
    Rt = rodrigues([0, 0, math.pi / 2])
    P = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [2, 1, 3]]
    Q = [_apply(Rt, [0, 0, 0], p) for p in P]
    R, t, rmsd = kabsch(P, Q)
    assert rmsd < 1e-6
    for a in range(3):
        for b in range(3):
            assert _close(R[a][b], Rt[a][b], 1e-5)


def test_length_mismatch_raises():
    with pytest.raises(ValueError):
        kabsch([[0, 0, 0]], [[0, 0, 0], [1, 1, 1]])


def test_wrong_dimension_raises():
    with pytest.raises(ValueError):
        kabsch([[0, 0]], [[0, 0]])


def test_empty_raises():
    with pytest.raises(ValueError):
        kabsch([], [])


def test_reflection_not_produced():
    # even for a point set whose best fit tempts a reflection, det stays +1
    rng = random.Random(473)
    P = [[rng.uniform(-1, 1) for _ in range(3)] for _ in range(6)]
    Q = [[-p[0], p[1], p[2]] for p in P]  # a reflection of P
    R, t, rmsd = kabsch(P, Q)
    assert _close(_det3(R), 1.0, 1e-6)  # proper rotation despite the reflected target
