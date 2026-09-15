"""Lanczos iteration: extremal eigenvalues of large symmetric matrices.

Full eigendecomposition (:func:`quantforge.pca.jacobi_eigen`) costs ``O(n^3)`` and needs the
dense matrix. When ``A`` is large and symmetric and you only want a few *extremal* eigenvalues
(the largest/smallest -- for spectral gaps, PCA of a huge covariance, graph connectivity), the
Lanczos iteration is far cheaper: it builds an orthonormal Krylov basis and projects ``A`` onto
it as a small ``m x m`` tridiagonal matrix ``T`` whose eigenvalues (the *Ritz values*) converge
to the extremal eigenvalues of ``A`` -- typically after ``m << n`` steps, using ``A`` only
through matrix-vector products.

This implements Lanczos with full reorthogonalization (stable for the moderate ``m`` used here)
and returns the Ritz values. ``A`` may be given as a matrix or a matrix-vector product callable.
Pure standard library on top of :func:`quantforge.pca.jacobi_eigen`.
"""

import math

from .pca import jacobi_eigen


def _as_matvec(A):
    if callable(A):
        return A
    return lambda v: [sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A))]


def _dot(a, b):
    return sum(a[i] * b[i] for i in range(len(a)))


def lanczos(A, n, m, v0=None):
    """Lanczos tridiagonalization of a symmetric operator ``A`` (dimension ``n``, ``m`` steps).

    ``A`` is an ``n x n`` symmetric matrix or a callable computing ``A @ v``. Returns
    ``(alpha, beta)`` -- the diagonal (length ``m``) and off-diagonal (length ``m-1``) of the
    tridiagonal projection ``T`` -- using full reorthogonalization against the stored basis.
    """
    matvec = _as_matvec(A)
    if v0 is None:
        v0 = [1.0 / math.sqrt(n)] * n
    else:
        nrm = math.sqrt(_dot(v0, v0))
        v0 = [x / nrm for x in v0]

    V = [v0]
    alpha = []
    beta = []
    w = matvec(v0)
    a = _dot(w, v0)
    alpha.append(a)
    w = [w[i] - a * v0[i] for i in range(n)]
    for j in range(1, m):
        # full reorthogonalization against all previous basis vectors
        for u in V:
            c = _dot(w, u)
            w = [w[i] - c * u[i] for i in range(n)]
        b = math.sqrt(_dot(w, w))
        if b < 1e-12:
            break                       # invariant subspace found
        beta.append(b)
        vj = [w[i] / b for i in range(n)]
        V.append(vj)
        w = matvec(vj)
        a = _dot(w, vj)
        alpha.append(a)
        w = [w[i] - a * vj[i] - b * V[-2][i] for i in range(n)]
    return alpha, beta


def _tridiag_eigenvalues(alpha, beta):
    # build the dense tridiagonal and use the symmetric eigensolver (m is small)
    m = len(alpha)
    T = [[0.0] * m for _ in range(m)]
    for i in range(m):
        T[i][i] = alpha[i]
    for i in range(m - 1):
        T[i][i + 1] = beta[i]
        T[i + 1][i] = beta[i]
    vals, _ = jacobi_eigen(T)
    return sorted(vals)


def lanczos_eigenvalues(A, n, m, v0=None):
    """Ritz eigenvalue estimates of symmetric ``A`` from ``m`` Lanczos steps (ascending).

    The extreme Ritz values (smallest and largest) converge fastest to the extremal eigenvalues
    of ``A``. Returns the sorted list of Ritz values (length up to ``m``).
    """
    alpha, beta = lanczos(A, n, m, v0)
    return _tridiag_eigenvalues(alpha, beta)
