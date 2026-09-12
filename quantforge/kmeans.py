"""K-means clustering (Lloyd's algorithm with k-means++ initialization).

Partitions ``n`` points into ``k`` clusters by alternating assignment (each point
to its nearest centroid) and update (each centroid to its cluster mean) until the
assignments stop changing. Centroids are seeded by k-means++, which spreads the
initial centers by sampling each with probability proportional to its squared
distance from the nearest chosen center -- far more reliable than random seeding.
A deterministic linear-congruential stream keeps runs reproducible per seed. Pure
standard library.
"""

import math


def _lcg(seed):
    state = seed & 0x7FFFFFFF

    def _next():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return state / 0x80000000

    return _next


def _dist2(a, b):
    return sum((a[i] - b[i]) ** 2 for i in range(len(a)))


def _kpp_init(X, k, rand):
    """k-means++ seeding: return k initial centroids."""
    n = len(X)
    first = int(rand() * n)
    centers = [list(X[first])]
    d2 = [_dist2(X[i], centers[0]) for i in range(n)]
    for _ in range(1, k):
        total = sum(d2)
        if total <= 0.0:
            # All remaining points coincide with a center; pick any.
            centers.append(list(X[int(rand() * n)]))
            continue
        target = rand() * total
        cum = 0.0
        chosen = n - 1
        for i in range(n):
            cum += d2[i]
            if cum >= target:
                chosen = i
                break
        centers.append(list(X[chosen]))
        for i in range(n):
            nd = _dist2(X[i], centers[-1])
            if nd < d2[i]:
                d2[i] = nd
    return centers


def kmeans(X, k, max_iter=100, seed=1234567):
    """Cluster ``X`` into ``k`` groups with Lloyd's algorithm (k-means++ seed).

    Parameters
    ----------
    X : list[list[float]]
        ``n`` points in ``d`` dimensions.
    k : int
        Number of clusters (>= 1, <= n).
    max_iter : int
        Maximum Lloyd iterations.
    seed : int
        Seed for the deterministic k-means++ sampler.

    Returns
    -------
    dict
        ``labels`` (cluster index per point), ``centroids``, ``inertia`` (total
        within-cluster squared distance), and ``iterations``.
    """
    n = len(X)
    if n == 0:
        raise ValueError("need at least one point")
    if k < 1 or k > n:
        raise ValueError("k must satisfy 1 <= k <= n")
    d = len(X[0])
    rand = _lcg(seed)
    centers = _kpp_init(X, k, rand)

    labels = [0] * n
    it = 0
    for it in range(1, max_iter + 1):
        changed = False
        for i in range(n):
            best, best_d = 0, float("inf")
            for c in range(k):
                dd = _dist2(X[i], centers[c])
                if dd < best_d:
                    best_d, best = dd, c
            if labels[i] != best:
                changed = True
            labels[i] = best
        # Update centroids to cluster means (keep empty clusters put).
        sums = [[0.0] * d for _ in range(k)]
        counts = [0] * k
        for i in range(n):
            c = labels[i]
            counts[c] += 1
            for j in range(d):
                sums[c][j] += X[i][j]
        for c in range(k):
            if counts[c] > 0:
                centers[c] = [sums[c][j] / counts[c] for j in range(d)]
        if not changed:
            break

    inertia = sum(_dist2(X[i], centers[labels[i]]) for i in range(n))
    return {
        "labels": labels,
        "centroids": centers,
        "inertia": inertia,
        "iterations": it,
    }
