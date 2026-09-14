"""Unit quaternions for 3-D rotation."""

import math
import random

import pytest

from quantforge import (
    quat_multiply,
    quat_normalize,
    axis_angle_to_quat,
    quat_to_axis_angle,
    rotate_vector,
    slerp,
)


def _close(a, b, t=1e-9):
    return all(abs(a[i] - b[i]) < t for i in range(len(a)))


def test_known_rotations():
    q = axis_angle_to_quat([0, 0, 1], math.pi / 2)      # 90 deg about z
    assert _close(rotate_vector(q, (1, 0, 0)), (0, 1, 0))
    q = axis_angle_to_quat([1, 0, 0], math.pi)          # 180 deg about x
    assert _close(rotate_vector(q, (0, 1, 0)), (0, -1, 0))


def test_rotation_preserves_length():
    rng = random.Random(1)
    for _ in range(500):
        axis = [rng.uniform(-1, 1) for _ in range(3)]
        if math.hypot(*axis) < 1e-6:
            continue
        q = axis_angle_to_quat(axis, rng.uniform(0, math.pi))
        v = (rng.uniform(-5, 5), rng.uniform(-5, 5), rng.uniform(-5, 5))
        rv = rotate_vector(q, v)
        assert abs(math.hypot(*rv) - math.hypot(*v)) < 1e-9


def test_axis_angle_roundtrip():
    rng = random.Random(2)
    for _ in range(500):
        axis = [rng.uniform(-1, 1) for _ in range(3)]
        n = math.hypot(*axis)
        if n < 1e-6:
            continue
        axis = [a / n for a in axis]
        angle = rng.uniform(0.01, math.pi - 0.01)
        ax2, a2 = quat_to_axis_angle(axis_angle_to_quat(axis, angle))
        assert abs(a2 - angle) < 1e-9
        assert _close(ax2, axis, 1e-7) or _close(ax2, [-a for a in axis], 1e-7)


def test_composition():
    q1 = axis_angle_to_quat([0, 0, 1], 0.5)
    q2 = axis_angle_to_quat([1, 0, 0], 0.3)
    v = (1, 2, 3)
    seq = rotate_vector(q2, rotate_vector(q1, v))
    comp = rotate_vector(quat_multiply(q2, q1), v)
    assert _close(seq, comp, 1e-9)


def test_slerp():
    a = axis_angle_to_quat([0, 0, 1], 0)
    b = axis_angle_to_quat([0, 0, 1], math.pi / 2)
    assert _close(slerp(a, b, 0), quat_normalize(a))
    mid = slerp(a, b, 0.5)                                # 45 deg rotation
    assert _close(rotate_vector(mid, (1, 0, 0)),
                  (math.cos(math.pi / 4), math.sin(math.pi / 4), 0), 1e-9)
    rng = random.Random(3)
    for _ in range(200):
        s = slerp(a, b, rng.random())
        assert abs(math.sqrt(sum(c * c for c in s)) - 1) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        quat_normalize((0, 0, 0, 0))
    with pytest.raises(ValueError):
        axis_angle_to_quat([0, 0, 0], 1)
