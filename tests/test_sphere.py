"""Tests for unit-sphere geometry, cross-checked against angle formulas."""

import math
import random

import pytest

from quantforge.sphere import (
    angular_distance,
    slerp_vectors,
    spherical_centroid,
    spherical_resultant_length,
)


def _close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def _norm(v):
    return math.sqrt(sum(c * c for c in v))


def _dot(a, b):
    return sum(a[i] * b[i] for i in range(3))


def test_fuzz_angular_vs_acos():
    rng = random.Random(481)
    for _ in range(5000):
        u = [rng.uniform(-5, 5) for _ in range(3)]
        v = [rng.uniform(-5, 5) for _ in range(3)]
        if _norm(u) < 1e-6 or _norm(v) < 1e-6:
            continue
        un = [c / _norm(u) for c in u]
        vn = [c / _norm(v) for c in v]
        d = max(-1, min(1, _dot(un, vn)))
        assert _close(angular_distance(u, v), math.acos(d), 1e-7)


def test_fuzz_slerp_endpoints_and_speed():
    rng = random.Random(482)
    for _ in range(5000):
        u = [rng.uniform(-5, 5) for _ in range(3)]
        v = [rng.uniform(-5, 5) for _ in range(3)]
        if _norm(u) < 1e-6 or _norm(v) < 1e-6:
            continue
        un = [c / _norm(u) for c in u]
        vn = [c / _norm(v) for c in v]
        assert all(_close(a, b, 1e-7) for a, b in zip(slerp_vectors(u, v, 0.0), un))
        assert all(_close(a, b, 1e-7) for a, b in zip(slerp_vectors(u, v, 1.0), vn))
        t = rng.random()
        st = slerp_vectors(u, v, t)
        assert _close(_norm(st), 1.0, 1e-9)
        ang = angular_distance(u, v)
        if 1e-6 < ang < math.pi - 1e-6:
            assert _close(angular_distance(un, st), t * ang, 1e-6)


def test_explicit_angles():
    assert _close(angular_distance([1, 0, 0], [0, 1, 0]), math.pi / 2)
    assert _close(angular_distance([1, 0, 0], [-1, 0, 0]), math.pi)
    assert _close(angular_distance([1, 0, 0], [1, 0, 0]), 0.0)


def test_slerp_midpoint():
    m = slerp_vectors([1, 0, 0], [0, 1, 0], 0.5)
    assert _close(m[0], 1 / math.sqrt(2), 1e-9)
    assert _close(m[1], 1 / math.sqrt(2), 1e-9)


def test_slerp_nearly_parallel():
    u = [1, 0, 0]
    v = [1, 1e-12, 0]
    s = slerp_vectors(u, v, 0.5)
    assert _close(_norm(s), 1.0)


def test_centroid_identical_directions():
    c = spherical_centroid([[1, 0, 0], [2, 0, 0], [5, 0, 0]])
    assert _close(c[0], 1.0)
    assert _close(c[1], 0.0)
    assert _close(c[2], 0.0)


def test_centroid_is_unit():
    c = spherical_centroid([[1, 1, 0], [0, 1, 1], [1, 0, 1]])
    assert _close(_norm(c), 1.0)


def test_resultant_bounds():
    assert _close(spherical_resultant_length([[1, 0, 0], [3, 0, 0]]), 1.0)  # coincident
    assert _close(spherical_resultant_length([[1, 0, 0], [-1, 0, 0]]), 0.0)  # opposite
    R = spherical_resultant_length([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert 0 < R < 1


def test_zero_vector_raises():
    with pytest.raises(ValueError):
        angular_distance([0, 0, 0], [1, 0, 0])


def test_empty_raises():
    with pytest.raises(ValueError):
        spherical_centroid([])
    with pytest.raises(ValueError):
        spherical_resultant_length([])


def test_cancelling_centroid_raises():
    with pytest.raises(ValueError):
        spherical_centroid([[1, 0, 0], [-1, 0, 0]])
