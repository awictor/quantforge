"""DBSCAN density-based clustering.

Unlike k-means, DBSCAN needs no preset cluster count, finds arbitrarily-shaped clusters,
and labels low-density points as *noise*. A point is a *core* point if at least
``min_samples`` points lie within radius ``eps``; clusters grow by connecting core
points and their neighbours, and anything not reachable from a core point is an outlier
(label ``-1``). Ideal for spatial data, anomaly detection, and any clustering where the
number of groups is unknown. Pure standard library; ``O(n^2)`` neighbour search.
"""

import math


def _dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(len(a))))


def dbscan(X, eps, min_samples):
    """Cluster points ``X`` by DBSCAN; return a label per point.

    ``eps`` is the neighbourhood radius, ``min_samples`` the core-point threshold
    (counting the point itself). Labels are ``0, 1, ...`` for clusters and ``-1`` for
    noise. No cluster count is required; clusters may be non-convex.
    """
    n = len(X)
    if n == 0:
        return []
    if eps <= 0:
        raise ValueError("eps must be positive")
    if min_samples < 1:
        raise ValueError("min_samples must be >= 1")

    # Precompute neighbourhoods.
    neighbors = [[j for j in range(n) if _dist(X[i], X[j]) <= eps] for i in range(n)]

    labels = [None] * n            # None = unvisited
    cluster = -1
    for i in range(n):
        if labels[i] is not None:
            continue
        if len(neighbors[i]) < min_samples:
            labels[i] = -1         # provisional noise (may be claimed as border later)
            continue
        cluster += 1
        labels[i] = cluster
        # Expand via a queue of density-reachable points.
        seeds = list(neighbors[i])
        k = 0
        while k < len(seeds):
            q = seeds[k]
            k += 1
            if labels[q] == -1:
                labels[q] = cluster      # border point
            if labels[q] is not None:
                continue
            labels[q] = cluster
            if len(neighbors[q]) >= min_samples:
                seeds.extend(neighbors[q])
    return labels
