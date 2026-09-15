"""Classical multidimensional scaling (MDS) -- embed points from a distance matrix.

Given only the pairwise distances between objects (not their coordinates), classical MDS finds a
low-dimensional point configuration whose Euclidean distances best match them. It is the metric
foundation of dimensionality reduction and the way to visualize similarity/dissimilarity data
(psychometrics, genetics, network layout). The recipe: square the distances, *double-center* the
matrix (``B = -1/2 J D2 J`` with ``J = I - 11^T/n``) so it becomes a Gram matrix of centred
coordinates, then take the top-``k`` eigenvectors scaled by the square roots of their
eigenvalues -- exactly the coordinates.

The embedding is determined only up to rotation, reflection, and translation, so it reproduces
the input distances (not the original coordinates). Uses :func:`quantforge.pca.jacobi_eigen`
(which returns eigenvectors as rows). Pure standard library.
"""

import math

from .pca import jacobi_eigen


def classical_mds(distances, k=2):
    """Classical (Torgerson) MDS: embed ``n`` objects in ``k`` dimensions from a distance matrix.

    ``distances`` is an ``n x n`` symmetric matrix of pairwise distances (zero diagonal). Returns
    a dict with ``coords`` (``n x k`` embedded points) and ``eigenvalues`` (the top ``k``, whose
    magnitudes measure how much each dimension explains). The configuration reproduces the input
    distances up to rotation/reflection/translation.
    """
    n = len(distances)
    if any(len(row) != n for row in distances):
        raise ValueError("distance matrix must be square")
    if not (1 <= k <= n):
        raise ValueError("require 1 <= k <= n")
    # squared distances
    D2 = [[distances[i][j] ** 2 for j in range(n)] for i in range(n)]
    # double-centering: B = -1/2 J D2 J  ==  -1/2 (D2 - row_mean - col_mean + grand_mean)
    row_mean = [sum(D2[i]) / n for i in range(n)]
    grand = sum(row_mean) / n
    col_mean = [sum(D2[i][j] for i in range(n)) / n for j in range(n)]
    B = [[-0.5 * (D2[i][j] - row_mean[i] - col_mean[j] + grand) for j in range(n)]
         for i in range(n)]

    vals, vecs = jacobi_eigen(B)              # eigenvectors are ROWS: vecs[t] is the t-th
    order = sorted(range(n), key=lambda t: -vals[t])[:k]   # largest eigenvalues
    coords = []
    for i in range(n):
        row = []
        for t in order:
            lam = vals[t]
            scale = math.sqrt(lam) if lam > 0 else 0.0
            row.append(vecs[t][i] * scale)
        coords.append(row)
    return {"coords": coords, "eigenvalues": [vals[t] for t in order]}


def embedded_distances(coords):
    """Pairwise Euclidean distance matrix of an ``n x k`` coordinate list (for validation)."""
    n = len(coords)
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = math.sqrt(sum((coords[i][t] - coords[j][t]) ** 2 for t in range(len(coords[i]))))
            D[i][j] = D[j][i] = d
    return D
