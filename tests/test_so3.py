"""Tests for the SO(3) exp/log map, cross-checked against the quaternion rotation path."""

import math
import random

import pytest

from quantforge.so3 import rodrigues, so3_log, hat, unhat
from quantforge.rotation import quat_to_matrix
from quantforge.quaternion import axis_angle_to_quat, quat_normalize


def _close(a, b, tol=1e-8):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def _mclose(A, B, tol=1e-8):
    return all(_close(A[i][j], B[i][j], tol) for i in range(3) for j in range(3))


def _vclose(a, b, tol=1e-7):
    return all(_close(x, y, tol) for x, y in zip(a, b))


def test_hat_unhat_inverse_and_cross():
    rng = random.Random(441)
    for _ in range(2000):
        v = [rng.uniform(-3, 3) for _ in range(3)]
        assert _vclose(unhat(hat(v)), v)
        w = [rng.uniform(-3, 3) for _ in range(3)]
        hw = [sum(hat(v)[i][j] * w[j] for j in range(3)) for i in range(3)]
        cross = [v[1] * w[2] - v[2] * w[1], v[2] * w[0] - v[0] * w[2], v[0] * w[1] - v[1] * w[0]]
        assert _vclose(hw, cross)


def test_fuzz_rodrigues_matches_quaternion_and_log_roundtrip():
    rng = random.Random(442)
    for _ in range(5000):
        axis = [rng.uniform(-1, 1) for _ in range(3)]
        n = math.sqrt(sum(a * a for a in axis))
        if n < 1e-9:
            continue
        axis = [a / n for a in axis]
        theta = rng.uniform(-math.pi + 0.01, math.pi - 0.01)
        omega = [theta * a for a in axis]
        R = rodrigues(omega)
        for i in range(3):
            for j in range(3):
                d = sum(R[i][k] * R[j][k] for k in range(3))
                assert _close(d, 1.0 if i == j else 0.0, 1e-7)
        q = quat_normalize(axis_angle_to_quat(axis, theta))
        assert _mclose(R, quat_to_matrix(q), 1e-7)
        rec = so3_log(R)
        assert _mclose(rodrigues(rec), R, 1e-7)


def test_90_about_z():
    R = rodrigues([0, 0, math.pi / 2])
    got = [sum(R[i][j] * [1, 0, 0][j] for j in range(3)) for i in range(3)]
    assert _vclose(got, [0, 1, 0])
    assert _vclose(so3_log(R), [0, 0, math.pi / 2])


def test_identity():
    assert _mclose(rodrigues([0, 0, 0]), [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert _vclose(so3_log([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), (0, 0, 0))


def test_pi_rotation():
    R = rodrigues([math.pi, 0, 0])
    lg = so3_log(R)
    assert _mclose(rodrigues(lg), R, 1e-6)


def test_small_angle():
    omega = [1e-10, 0, 0]
    R = rodrigues(omega)
    assert _mclose(R, [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 1e-8)


def test_log_of_small_rotation():
    omega = [1e-4, 2e-4, -1e-4]
    R = rodrigues(omega)
    assert _vclose(so3_log(R), omega, 1e-6)


def test_composition_is_rotation():
    # exp(a) exp(b) is still a rotation whose log round-trips
    Ra = rodrigues([0.3, -0.2, 0.5])
    Rb = rodrigues([0.1, 0.4, -0.3])
    Rc = [[sum(Ra[i][k] * Rb[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    assert _mclose(rodrigues(so3_log(Rc)), Rc, 1e-7)


def test_hat_is_skew_symmetric():
    K = hat([1, 2, 3])
    for i in range(3):
        for j in range(3):
            assert _close(K[i][j], -K[j][i])
