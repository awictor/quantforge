"""Silhouette score for clustering quality.

For each point, the silhouette compares how close it sits to its own cluster
versus the nearest other cluster:

    s(i) = (b(i) - a(i)) / max(a(i), b(i)),

where ``a(i)`` is the mean distance to the other points in its cluster and
``b(i)`` is the smallest mean distance to any other cluster. ``s`` ranges over
``[-1, 1]``: near 1 = well inside its cluster, near 0 = on a boundary, negative =
probably misassigned. The overall silhouette score is the mean of ``s(i)``, a
label-free way to compare clusterings (e.g. to choose ``k``). Pure standard
library.
"""

import math


def _dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(len(a))))


def silhouette_samples(X, labels):
    """Per-point silhouette values ``s(i)`` in ``[-1, 1]``.

    A point alone in its cluster gets ``s = 0`` by convention. Requires at least
    two distinct cluster labels.
    """
    n = len(X)
    if n != len(labels):
        raise ValueError("X and labels must have the same length")
    if n == 0:
        raise ValueError("need at least one point")
    clusters = sorted(set(labels))
    if len(clusters) < 2:
        raise ValueError("silhouette needs at least 2 clusters")
    members = {c: [i for i in range(n) if labels[i] == c] for c in clusters}

    out = []
    for i in range(n):
        own = labels[i]
        same = members[own]
        if len(same) <= 1:
            out.append(0.0)
            continue
        a = sum(_dist(X[i], X[j]) for j in same if j != i) / (len(same) - 1)
        b = min(
            sum(_dist(X[i], X[j]) for j in members[c]) / len(members[c])
            for c in clusters if c != own
        )
        denom = max(a, b)
        out.append((b - a) / denom if denom > 0 else 0.0)
    return out


def silhouette_score(X, labels):
    """Mean silhouette over all points -- an overall clustering-quality score.

    Near 1 = dense, well-separated clusters; near 0 = overlapping; negative =
    mostly misassigned. Use it to compare label sets or pick ``k``.
    """
    s = silhouette_samples(X, labels)
    return sum(s) / len(s)
