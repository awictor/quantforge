"""Spectral clustering: cluster by the eigenvectors of the graph Laplacian.

k-means (:func:`quantforge.kmeans.kmeans`) assumes convex, roughly spherical clusters and fails
on curved or intertwined structure (nested rings, moons). Spectral clustering sidesteps that: it
builds a similarity graph over the points, forms its (normalized) Laplacian, embeds the points in
the space of the Laplacian's smallest eigenvectors -- where the manifold structure unfolds into
tight blobs -- and runs k-means there. It is the standard method for non-convex clusters and the
graph-partitioning view of clustering (normalized cut).

This uses an RBF (Gaussian) affinity, the symmetric normalized Laplacian
``L = I - D^{-1/2} W D^{-1/2}``, its smallest ``k`` eigenvectors (from
:func:`quantforge.pca.jacobi_eigen`), row-normalized (Ng-Jordan-Weiss), and
:func:`quantforge.kmeans.kmeans`. Pure standard library.
"""

import math

from .pca import jacobi_eigen
from .kmeans import kmeans


def rbf_affinity(X, gamma=1.0):
    """RBF (Gaussian) affinity matrix ``W[i][j] = exp(-gamma ||x_i - x_j||^2)`` (diagonal 1)."""
    n = len(X)
    d = len(X[0])
    W = [[0.0] * n for _ in range(n)]
    for i in range(n):
        W[i][i] = 1.0
        for j in range(i + 1, n):
            dist2 = sum((X[i][t] - X[j][t]) ** 2 for t in range(d))
            w = math.exp(-gamma * dist2)
            W[i][j] = W[j][i] = w
    return W


def spectral_clustering(X, k, gamma=1.0, seed=1234567):
    """Cluster points ``X`` into ``k`` groups by spectral clustering (Ng-Jordan-Weiss).

    Builds an RBF affinity with bandwidth ``gamma``, forms the symmetric normalized Laplacian,
    embeds the points in its ``k`` smallest eigenvectors (row-normalized), and k-means-clusters
    the embedding. Returns a dict with ``labels`` (cluster index per point) and ``embedding``.
    """
    n = len(X)
    if not (1 <= k <= n):
        raise ValueError("require 1 <= k <= n points")
    W = rbf_affinity(X, gamma)
    deg = [sum(W[i]) for i in range(n)]
    dinv = [1.0 / math.sqrt(deg[i]) if deg[i] > 0 else 0.0 for i in range(n)]
    # symmetric normalized Laplacian L = I - D^-1/2 W D^-1/2
    L = [[(1.0 if i == j else 0.0) - dinv[i] * W[i][j] * dinv[j] for j in range(n)]
         for i in range(n)]
    vals, vecs = jacobi_eigen(L)               # jacobi_eigen returns eigenvectors as ROWS
    # take the k eigenvectors of the smallest eigenvalues
    order = sorted(range(n), key=lambda t: vals[t])[:k]
    # embedding rows, then row-normalize (NJW); vecs[t] is the t-th eigenvector
    emb = []
    for i in range(n):
        row = [vecs[t][i] for t in order]
        nrm = math.sqrt(sum(v * v for v in row)) or 1.0
        emb.append([v / nrm for v in row])
    res = kmeans(emb, k, seed=seed)
    return {"labels": res["labels"], "embedding": emb}
