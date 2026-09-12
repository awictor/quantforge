"""Principal component analysis via the Jacobi eigenvalue algorithm.

Symmetric-matrix eigendecomposition (Jacobi rotations) and PCA of a covariance
matrix -- the standard yield-curve factor decomposition into level, slope, and
curvature. Returns eigenvalues (component variances) sorted descending with their
orthonormal eigenvectors (loadings), plus the proportion of variance explained.
Pure standard library.
"""

import math


def jacobi_eigen(matrix, tol=1e-12, max_sweeps=100):
    """Eigenvalues and eigenvectors of a symmetric matrix (Jacobi rotations).

    Returns ``(eigenvalues, eigenvectors)`` where ``eigenvectors[i]`` is the
    orthonormal eigenvector for ``eigenvalues[i]``, sorted by descending
    eigenvalue. Requires a symmetric input; iteratively zeroes off-diagonal
    entries with plane rotations.
    """
    n = len(matrix)
    for i in range(n):
        if len(matrix[i]) != n:
            raise ValueError("matrix must be square")
        for j in range(n):
            if abs(matrix[i][j] - matrix[j][i]) > 1e-9:
                raise ValueError("matrix must be symmetric")
    a = [row[:] for row in matrix]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(max_sweeps):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(i + 1, n))
        if off < tol:
            break
        for p in range(n):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-300:
                    continue
                theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
                t = (1.0 if theta >= 0 else -1.0) / \
                    (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq
    eigenvalues = [a[i][i] for i in range(n)]
    eigenvectors = [[v[r][c] for r in range(n)] for c in range(n)]
    order = sorted(range(n), key=lambda i: eigenvalues[i], reverse=True)
    eigenvalues = [eigenvalues[i] for i in order]
    eigenvectors = [eigenvectors[i] for i in order]
    return eigenvalues, eigenvectors


def pca(covariance):
    """PCA of a covariance matrix: sorted variances, loadings, variance explained.

    Returns a dict with ``variances`` (eigenvalues, descending), ``loadings``
    (orthonormal eigenvectors), ``explained`` (each variance over the total), and
    ``cumulative_explained``. For a yield-curve covariance the first three
    components are the level, slope, and curvature factors.
    """
    variances, loadings = jacobi_eigen(covariance)
    total = sum(variances)
    if total <= 0:
        raise ValueError("covariance must have positive total variance")
    explained = [v / total for v in variances]
    cum = []
    acc = 0.0
    for e in explained:
        acc += e
        cum.append(acc)
    return {"variances": variances, "loadings": loadings,
            "explained": explained, "cumulative_explained": cum}


def reconstruct_covariance(variances, loadings, k=None):
    """Rebuild a covariance matrix from the top ``k`` principal components.

    ``sum_i variance_i * (loading_i outer loading_i)`` over the first ``k``
    components. With all components it reproduces the original covariance exactly
    (spectral decomposition); with ``k`` below the rank it is the best rank-``k``
    approximation.
    """
    n = len(loadings)
    m = len(loadings[0]) if n else 0
    k = n if k is None else min(k, n)
    out = [[0.0] * m for _ in range(m)]
    for i in range(k):
        lam = variances[i]
        vec = loadings[i]
        for r in range(m):
            for c in range(m):
                out[r][c] += lam * vec[r] * vec[c]
    return out


def pca_scenario(component_index, n_sigma, variances, loadings):
    """A stress scenario shocking one principal component by ``n_sigma`` std devs.

    Returns the vector move ``n_sigma * sqrt(variance_i) * loading_i`` -- a
    ``n_sigma``-standard-deviation move along principal component
    ``component_index``. For a yield curve, component 0 is a parallel (level)
    shift, 1 a slope twist, 2 a curvature bend.
    """
    n = len(loadings)
    if not (0 <= component_index < n):
        raise ValueError("component_index out of range")
    if variances[component_index] < 0:
        raise ValueError("variance must be non-negative")
    scale = n_sigma * math.sqrt(variances[component_index])
    return [scale * x for x in loadings[component_index]]


def project(data_row, loadings, k=None):
    """Project a data vector onto the first ``k`` principal components (scores).

    ``sum_j data_row_j loadings_i_j`` for each retained component ``i``. ``k``
    defaults to all components. The scores are the coordinates of the observation
    in the principal-component basis.
    """
    n = len(loadings)
    k = n if k is None else min(k, n)
    if any(len(load) != len(data_row) for load in loadings):
        raise ValueError("loadings and data dimensions must match")
    return [sum(data_row[j] * loadings[i][j] for j in range(len(data_row)))
            for i in range(k)]
