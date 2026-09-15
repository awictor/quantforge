"""Isomap: nonlinear dimensionality reduction by geodesic MDS.

Classical MDS (:mod:`quantforge.mds`) embeds points to preserve *straight-line* distances, which
flattens a curved manifold incorrectly -- points on opposite ends of a rolled-up sheet look close
in Euclidean distance but are far along the surface. Isomap (Tenenbaum, de Silva & Langford 2000)
fixes this: it builds a k-nearest-neighbour graph, approximates the *geodesic* (along-manifold)
distance between every pair as the shortest path through that graph, and then runs classical MDS
on those geodesic distances. The result unrolls manifolds like the Swiss roll that linear methods
(PCA, classical MDS) cannot.

Uses :func:`quantforge.graph4.floyd_warshall` for all-pairs geodesics and
:func:`quantforge.mds.classical_mds` for the final embedding. Pure standard library.
"""

import math

from .graph4 import floyd_warshall
from .mds import classical_mds


def _euclidean(a, b):
    return math.sqrt(sum((a[t] - b[t]) ** 2 for t in range(len(a))))


def knn_graph(X, n_neighbors):
    """Symmetric k-nearest-neighbour graph of points ``X`` with Euclidean edge weights.

    Returns ``{i: {j: dist}}``; the graph is symmetrized (an edge is kept if either endpoint has
    the other among its ``n_neighbors`` nearest). Used as the connectivity for geodesic distances.
    """
    n = len(X)
    if not (1 <= n_neighbors < n):
        raise ValueError("require 1 <= n_neighbors < n")
    g = {i: {} for i in range(n)}
    for i in range(n):
        dists = sorted(((_euclidean(X[i], X[j]), j) for j in range(n) if j != i))
        for d, j in dists[:n_neighbors]:
            g[i][j] = d
            g[j][i] = d       # symmetrize
    return g


def isomap(X, k=2, n_neighbors=5):
    """Isomap embedding of points ``X`` into ``k`` dimensions.

    Builds the ``n_neighbors``-NN graph, computes all-pairs geodesic (shortest-path) distances,
    and runs classical MDS on them. Returns a dict with ``coords`` (``n x k``) and ``eigenvalues``.
    Raises if the neighbourhood graph is disconnected (some geodesic distance is infinite).
    """
    n = len(X)
    g = knn_graph(X, n_neighbors)
    dist = floyd_warshall(g)
    # dense geodesic distance matrix
    D = [[dist[i][j] for j in range(n)] for i in range(n)]
    if any(D[i][j] == float("inf") for i in range(n) for j in range(n)):
        raise ValueError("neighbourhood graph is disconnected; increase n_neighbors")
    return classical_mds(D, k)
