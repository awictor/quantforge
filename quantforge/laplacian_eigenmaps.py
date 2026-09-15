"""Laplacian eigenmaps: nonlinear dimensionality reduction by the graph Laplacian.

Isomap (:mod:`quantforge.isomap`) preserves *global* geodesic distances; Laplacian eigenmaps
(Belkin & Niyogi 2003) instead preserve *local* neighbourhood structure -- points close on the
manifold map to nearby embedding coordinates. It builds a k-nearest-neighbour graph with
heat-kernel weights, forms the graph Laplacian, and uses its smallest *non-trivial* eigenvectors
as the embedding. Minimizing ``sum_ij W_ij ||y_i - y_j||^2`` (which the Laplacian's Rayleigh
quotient does) keeps neighbours together, so the method excels at unfolding manifolds while
respecting local geometry -- the embedding step shared with spectral clustering.

The trivial constant eigenvector (eigenvalue 0) is dropped; the next ``k`` give the coordinates.
Uses :func:`quantforge.pca.jacobi_eigen` (which returns eigenvectors as ROWS). Pure standard
library.
"""

import math

from .pca import jacobi_eigen


def _euclidean2(a, b):
    return sum((a[t] - b[t]) ** 2 for t in range(len(a)))


def heat_knn_affinity(X, n_neighbors, t=1.0):
    """Symmetric k-NN affinity with heat-kernel weights ``exp(-||x_i - x_j||^2 / t)``.

    An edge ``(i, j)`` is kept if either point is among the other's ``n_neighbors`` nearest;
    its weight is the heat kernel of the squared distance. Returns an ``n x n`` matrix.
    """
    n = len(X)
    if not (1 <= n_neighbors < n):
        raise ValueError("require 1 <= n_neighbors < n")
    W = [[0.0] * n for _ in range(n)]
    for i in range(n):
        d = sorted(((_euclidean2(X[i], X[j]), j) for j in range(n) if j != i))
        for dist2, j in d[:n_neighbors]:
            w = math.exp(-dist2 / t)
            W[i][j] = w
            W[j][i] = w        # symmetrize
    return W


def laplacian_eigenmaps(X, k=2, n_neighbors=5, t=1.0):
    """Laplacian-eigenmap embedding of points ``X`` into ``k`` dimensions.

    Builds the heat-kernel k-NN affinity, forms the symmetric normalized Laplacian, and returns
    the ``k`` smallest *non-trivial* eigenvectors (dropping the constant eigenvector) as the
    embedding. Returns a dict with ``coords`` (``n x k``) and ``eigenvalues`` (the used ones).
    """
    n = len(X)
    if not (1 <= k <= n - 1):
        raise ValueError("require 1 <= k <= n - 1")
    W = heat_knn_affinity(X, n_neighbors, t)
    deg = [sum(W[i]) for i in range(n)]
    dinv = [1.0 / math.sqrt(deg[i]) if deg[i] > 0 else 0.0 for i in range(n)]
    # symmetric normalized Laplacian
    L = [[(1.0 if i == j else 0.0) - dinv[i] * W[i][j] * dinv[j] for j in range(n)]
         for i in range(n)]
    vals, vecs = jacobi_eigen(L)              # eigenvectors are ROWS
    order = sorted(range(n), key=lambda idx: vals[idx])
    # drop the first (trivial, eigenvalue ~0) eigenvector; take the next k
    chosen = order[1:1 + k]
    coords = [[vecs[t_][i] for t_ in chosen] for i in range(n)]
    return {"coords": coords, "eigenvalues": [vals[t_] for t_ in chosen]}
