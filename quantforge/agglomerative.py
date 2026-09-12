"""Agglomerative hierarchical clustering.

Builds a cluster tree bottom-up: start with every point its own cluster, then
repeatedly merge the two closest clusters until one remains. Inter-cluster
distance is the linkage:

  * ``single``   -- nearest pair (chains, finds elongated shapes).
  * ``complete`` -- farthest pair (compact, equal-diameter clusters).
  * ``average``  -- mean pairwise distance (a balance of the two).

``linkage`` returns the merge history with monotone (for single/complete/average)
merge distances; ``fcluster`` cuts the tree at a distance threshold or into a
target number of clusters. Pure standard library.
"""

import math


def _dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(len(a))))


def linkage(X, method="average"):
    """Agglomerative linkage over points ``X``.

    Returns a list of ``n-1`` merges, each ``(cluster_a, cluster_b, distance,
    size)``. Cluster ids ``0..n-1`` are the singletons; merge ``m`` creates the new
    id ``n + m``. Merge distances are non-decreasing for single/complete/average
    linkage.
    """
    if method not in ("single", "complete", "average"):
        raise ValueError("method must be single, complete, or average")
    n = len(X)
    if n < 2:
        raise ValueError("need at least 2 points")

    # Active clusters: id -> list of member point indices.
    active = {i: [i] for i in range(n)}
    next_id = n
    merges = []

    def cluster_dist(a_members, b_members):
        ds = [_dist(X[i], X[j]) for i in a_members for j in b_members]
        if method == "single":
            return min(ds)
        if method == "complete":
            return max(ds)
        return sum(ds) / len(ds)          # average

    while len(active) > 1:
        ids = list(active)
        best = None
        best_d = float("inf")
        for a in range(len(ids)):
            for b in range(a + 1, len(ids)):
                d = cluster_dist(active[ids[a]], active[ids[b]])
                if d < best_d:
                    best_d = d
                    best = (ids[a], ids[b])
        ca, cb = best
        merged = active[ca] + active[cb]
        merges.append((ca, cb, best_d, len(merged)))
        del active[ca]
        del active[cb]
        active[next_id] = merged
        next_id += 1
    return merges


def fcluster(X, merges, n_clusters=None, distance_threshold=None):
    """Flatten a linkage into cluster labels.

    Provide either ``n_clusters`` (cut so that many clusters remain) or
    ``distance_threshold`` (merge only below that distance). Returns a label per
    original point, relabeled to ``0..k-1`` in order of first appearance.
    """
    n = len(X)
    if (n_clusters is None) == (distance_threshold is None):
        raise ValueError("give exactly one of n_clusters or distance_threshold")

    # Replay merges, stopping per the criterion.
    parent = {i: [i] for i in range(n)}
    next_id = n
    if n_clusters is not None:
        if not (1 <= n_clusters <= n):
            raise ValueError("n_clusters must be in [1, n]")
        stop_after = n - n_clusters
        for k, (ca, cb, d, _) in enumerate(merges):
            if k >= stop_after:
                break
            parent[next_id] = parent[ca] + parent[cb]
            del parent[ca]
            del parent[cb]
            next_id += 1
    else:
        for (ca, cb, d, _) in merges:
            if d > distance_threshold:
                break
            parent[next_id] = parent[ca] + parent[cb]
            del parent[ca]
            del parent[cb]
            next_id += 1

    labels = [0] * n
    for lab, members in enumerate(parent.values()):
        for i in members:
            labels[i] = lab
    # Relabel by first appearance for determinism.
    remap = {}
    out = []
    for lab in labels:
        if lab not in remap:
            remap[lab] = len(remap)
        out.append(remap[lab])
    return out
