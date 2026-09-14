"""Tests for 3-D rotation conversions, cross-checked against quaternion.rotate_vector."""

import math
import random

import pytest

from quantforge.rotation import (
    quat_to_matrix,
    matrix_to_quat,
    euler_to_quat,
    quat_to_euler,
    euler_to_matrix,
    matrix_to_euler,
)
from quantforge.quaternion import rotate_vector, quat_normalize, axis_angle_to_quat


def _close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def _qclose(a, b, tol=1e-8):
    return all(_close(x, y, tol) for x, y in zip(a, b)) or all(_close(x, -y, tol) for x, y in zip(a, b))


def _matvec(M, v):
    return [sum(M[i][j] * v[j] for j in range(3)) for i in range(3)]


def _random_quat(rng):
    axis = [rng.uniform(-1, 1) for _ in range(3)]
    if all(a == 0 for a in axis):
        axis = [1, 0, 0]
    return quat_normalize(axis_angle_to_quat(axis, rng.uniform(-math.pi, math.pi)))


def test_fuzz_matrix_orthogonal_and_roundtrip():
    rng = random.Random(431)
    for _ in range(5000):
        q = _random_quat(rng)
        M = quat_to_matrix(q)
        for i in range(3):
            for j in range(3):
                dot = sum(M[i][k] * M[j][k] for k in range(3))
                assert _close(dot, 1.0 if i == j else 0.0, 1e-8)
        assert _qclose(matrix_to_quat(M), q)


def test_fuzz_matrix_matches_rotate_vector():
    rng = random.Random(432)
    for _ in range(5000):
        q = _random_quat(rng)
        M = quat_to_matrix(q)
        v = [rng.uniform(-3, 3) for _ in range(3)]
        mv = _matvec(M, v)
        rv = rotate_vector(q, v)
        assert all(_close(a, b, 1e-7) for a, b in zip(mv, rv))


def test_fuzz_euler_roundtrip():
    rng = random.Random(433)
    for _ in range(5000):
        yaw = rng.uniform(-3, 3)
        pitch = rng.uniform(-1.4, 1.4)  # avoid the gimbal-lock pole
        roll = rng.uniform(-3, 3)
        y2, p2, r2 = quat_to_euler(euler_to_quat(yaw, pitch, roll))
        assert _close(yaw, y2, 1e-6)
        assert _close(pitch, p2, 1e-6)
        assert _close(roll, r2, 1e-6)
        y3, p3, r3 = matrix_to_euler(euler_to_matrix(yaw, pitch, roll))
        assert _close(yaw, y3, 1e-6)
        assert _close(pitch, p3, 1e-6)
        assert _close(roll, r3, 1e-6)


def test_identity():
    assert quat_to_matrix((1, 0, 0, 0)) == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    assert quat_to_euler((1, 0, 0, 0)) == (0.0, 0.0, 0.0)


def test_90_degree_z_rotation():
    q = (math.cos(math.pi / 4), 0, 0, math.sin(math.pi / 4))
    M = quat_to_matrix(q)
    assert all(_close(a, b, 1e-9) for a, b in zip(_matvec(M, [1, 0, 0]), [0, 1, 0]))


def test_yaw_90_via_euler():
    M = euler_to_matrix(math.pi / 2, 0, 0)
    assert all(_close(a, b, 1e-9) for a, b in zip(_matvec(M, [1, 0, 0]), [0, 1, 0]))


def test_matrix_to_quat_identity():
    assert _qclose(matrix_to_quat([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), (1, 0, 0, 0))


def test_euler_pitch_only():
    # pitch 90 deg maps x -> -z (Z-Y-X)
    M = euler_to_matrix(0, math.pi / 2, 0)
    got = _matvec(M, [1, 0, 0])
    assert all(_close(a, b, 1e-9) for a, b in zip(got, [0, 0, -1]))


def test_zero_quaternion_raises():
    with pytest.raises(ValueError):
        quat_to_matrix((0, 0, 0, 0))
    with pytest.raises(ValueError):
        quat_to_euler((0, 0, 0, 0))


def test_matrix_to_quat_180_degree():
    # 180 about x: diag(1,-1,-1)
    M = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
    q = matrix_to_quat(M)
    # should be (0, 1, 0, 0) up to sign
    assert _qclose(q, (0, 1, 0, 0))
