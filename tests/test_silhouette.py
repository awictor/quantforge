"""Silhouette clustering-quality score."""

import random

import pytest

from quantforge import silhouette_score, silhouette_samples
from quantforge import kmeans


def _blobs(seed=1):
    rng = random.Random(seed)
    X = []
    for cx, cy in [(0, 0), (10, 10), (0, 10)]:
        for _ in range(40):
            X.append([cx + rng.gauss(0, 0.4), cy + rng.gauss(0, 0.4)])
    return X


def test_correct_k_high_score():
    X = _blobs()
    labels = kmeans(X, 3, seed=42)["labels"]
    assert silhouette_score(X, labels) > 0.7


def test_correct_k_beats_wrong_k():
    X = _blobs()
    s3 = silhouette_score(X, kmeans(X, 3, seed=42)["labels"])
    s2 = silhouette_score(X, kmeans(X, 2, seed=42)["labels"])
    assert s3 > s2


def test_samples_in_range():
    X = _blobs()
    labels = kmeans(X, 3, seed=42)["labels"]
    assert all(-1.0 <= v <= 1.0 for v in silhouette_samples(X, labels))


def test_singleton_cluster_scores_zero():
    X = _blobs()
    labels = kmeans(X, 3, seed=42)["labels"]
    labels = list(labels)
    labels[0] = 99                                   # isolate point 0
    assert silhouette_samples(X, labels)[0] == 0.0


def test_tight_separated_near_one():
    X = [[0, 0], [0, 0.01], [100, 100], [100, 100.01]]
    assert silhouette_score(X, [0, 0, 1, 1]) > 0.99


def test_random_labels_low():
    X = _blobs()
    rng = random.Random(2)
    labels = [rng.randint(0, 2) for _ in range(len(X))]
    assert silhouette_score(X, labels) < 0.2


def test_validation():
    X = _blobs()
    with pytest.raises(ValueError):
        silhouette_score(X, [0] * len(X))            # one cluster only
    with pytest.raises(ValueError):
        silhouette_score(X, [0, 1])                  # length mismatch
