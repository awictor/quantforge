"""Vector algebra: dot, cross, norms, angles, projection, reflection."""

import math
import random

import pytest

from quantforge import dot, cross, norm, normalize, angle_between, vector_project, vector_reject, reflect


def _close(a, b, t=1e-9):
    return all(abs(a[i] - b[i]) < t for i in range(len(a)))


def test_cross_basis_and_anticommute():
    assert cross((1, 0, 0), (0, 1, 0)) == (0, 0, 1)
    rng = random.Random(1)
    for _ in range(300):
        a = tuple(rng.uniform(-5, 5) for _ in range(3))
        b = tuple(rng.uniform(-5, 5) for _ in range(3))
        axb = cross(a, b)
        assert _close(axb, tuple(-c for c in cross(b, a)))
        assert abs(dot(axb, a)) < 1e-9 and abs(dot(axb, b)) < 1e-9   # perpendicular


def test_angles():
    assert abs(math.degrees(angle_between((1, 0, 0), (0, 1, 0))) - 90) < 1e-9
    assert abs(math.degrees(angle_between((1, 2, 3), (2, 4, 6)))) < 1e-6
    assert abs(math.degrees(angle_between((1, 0), (-1, 0))) - 180) < 1e-9


def test_projection_rejection():
    rng = random.Random(2)
    for _ in range(300):
        a = tuple(rng.uniform(-5, 5) for _ in range(3))
        b = tuple(rng.uniform(-5, 5) for _ in range(3))
        if norm(b) < 1e-6:
            continue
        p, r = vector_project(a, b), vector_reject(a, b)
        assert _close(tuple(p[i] + r[i] for i in range(3)), a)
        assert abs(dot(r, b)) < 1e-9


def test_reflect():
    assert _close(reflect((1, 2, 3), (0, 0, 1)), (1, 2, -3))       # flip z
    rng = random.Random(3)
    for _ in range(200):
        a = tuple(rng.uniform(-5, 5) for _ in range(3))
        n = tuple(rng.uniform(-1, 1) for _ in range(3))
        if norm(n) < 1e-6:
            continue
        assert abs(norm(reflect(a, n)) - norm(a)) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        cross((1, 2), (3, 4))
    with pytest.raises(ValueError):
        normalize((0, 0, 0))
    with pytest.raises(ValueError):
        dot((1, 2), (1, 2, 3))
