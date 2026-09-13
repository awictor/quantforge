"""DBSCAN density clustering."""

import math
import random

import pytest

from quantforge import dbscan


def _blobs(seed=3):
    rng = random.Random(seed)
    b1 = [[rng.gauss(0, 0.3), rng.gauss(0, 0.3)] for _ in range(50)]
    b2 = [[rng.gauss(10, 0.3), rng.gauss(10, 0.3)] for _ in range(50)]
    return b1 + b2


def test_two_blobs():
    labels = dbscan(_blobs(), eps=1.0, min_samples=5)
    assert len(set(l for l in labels if l >= 0)) == 2
    assert labels.count(-1) == 0
    assert len(set(labels[:50])) == 1 and len(set(labels[50:])) == 1
    assert set(labels[:50]) != set(labels[50:])


def test_outlier_is_noise():
    rng = random.Random(3)
    b1 = [[rng.gauss(0, 0.3), rng.gauss(0, 0.3)] for _ in range(50)]
    labels = dbscan(b1 + [[100.0, 100.0]], eps=1.0, min_samples=5)
    assert labels[-1] == -1


def test_non_convex_rings():
    inner = [[math.cos(t), math.sin(t)] for t in [i * 0.2 for i in range(32)]]
    outer = [[5 * math.cos(t), 5 * math.sin(t)] for t in [i * 0.2 for i in range(32)]]
    labels = dbscan(inner + outer, eps=0.8, min_samples=2)
    assert len(set(l for l in labels if l >= 0)) == 2


def test_tiny_eps_all_noise():
    assert all(l == -1 for l in dbscan(_blobs(), eps=0.001, min_samples=5))


def test_huge_eps_one_cluster():
    labels = dbscan(_blobs(), eps=100, min_samples=5)
    assert len(set(l for l in labels if l >= 0)) == 1


def test_deterministic():
    X = _blobs()
    assert dbscan(X, eps=1.0, min_samples=5) == dbscan(X, eps=1.0, min_samples=5)


def test_validation():
    with pytest.raises(ValueError):
        dbscan([[0, 0]], eps=0, min_samples=1)
    with pytest.raises(ValueError):
        dbscan([[0, 0]], eps=1, min_samples=0)
