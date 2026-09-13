"""Singular value decomposition by the one-sided Jacobi method, and the pseudo-inverse.

Factors an ``m x n`` matrix (``m >= n``) as ``A = U S V'`` with ``U`` (``m x n``)
and ``V`` (``n x n``) having orthonormal columns and ``S`` the non-negative singular
values in descending order. The one-sided Jacobi scheme orthogonalizes the columns
of ``A`` by a sequence of plane rotations -- simple, and accurate even for the small
singular values. The Moore-Penrose pseudo-inverse ``A^+ = V S^+ U'`` follows
directly and gives the minimum-norm least-squares solution. Pure standard library.
"""

import math


def svd(A, tol=1e-14, max_sweeps=60):
    """One-sided Jacobi SVD of an ``m x n`` matrix (``m >= n``).

    Returns ``(U, s, V)`` where ``U`` is ``m x n`` with orthonormal columns, ``s`` is
    the length-``n`` list of singular values (descending, non-negative), and ``V`` is
    ``n x n`` orthogonal, such that ``A = U diag(s) V'``. Iteratively rotates column
    pairs until they are orthogonal.
    """
    m = len(A)
    if m == 0:
        raise ValueError("A must be non-empty")
    n = len(A[0])
    if any(len(row) != n for row in A):
        raise ValueError("A must be rectangular")
    if m < n:
        raise ValueError("need m >= n rows")

    # Work on a copy of the columns of A; V accumulates the rotations.
    U = [row[:] for row in A]
    V = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    for _ in range(max_sweeps):
        off = 0.0
        for p in range(n - 1):
            for q in range(p + 1, n):
                alpha = sum(U[i][p] * U[i][p] for i in range(m))
                beta = sum(U[i][q] * U[i][q] for i in range(m))
                gamma = sum(U[i][p] * U[i][q] for i in range(m))
                off += gamma * gamma
                if abs(gamma) < 1e-300:
                    continue
                zeta = (beta - alpha) / (2.0 * gamma)
                t = (1.0 if zeta >= 0 else -1.0) / (abs(zeta) + math.sqrt(1.0 + zeta * zeta))
                c = 1.0 / math.sqrt(1.0 + t * t)
                sn = c * t
                for i in range(m):
                    up, uq = U[i][p], U[i][q]
                    U[i][p] = c * up - sn * uq
                    U[i][q] = sn * up + c * uq
                for i in range(n):
                    vp, vq = V[i][p], V[i][q]
                    V[i][p] = c * vp - sn * vq
                    V[i][q] = sn * vp + c * vq
        if off < tol:
            break

    # Singular values are the column norms of the orthogonalized U; normalize U.
    s = [math.sqrt(sum(U[i][j] * U[i][j] for i in range(m))) for j in range(n)]
    for j in range(n):
        if s[j] > 1e-300:
            for i in range(m):
                U[i][j] /= s[j]
    # Sort by descending singular value.
    order = sorted(range(n), key=lambda j: s[j], reverse=True)
    s_sorted = [s[j] for j in order]
    U_sorted = [[U[i][order[j]] for j in range(n)] for i in range(m)]
    V_sorted = [[V[i][order[j]] for j in range(n)] for i in range(n)]
    return U_sorted, s_sorted, V_sorted


def pseudo_inverse(A, rcond=1e-12):
    """Moore-Penrose pseudo-inverse ``A^+`` via the SVD.

    ``A^+ = V S^+ U'`` with the reciprocals of the singular values above
    ``rcond * s_max`` (smaller ones treated as zero). For a full-rank tall ``A`` this
    is ``(A' A)^{-1} A'``; applied to ``b`` it gives the minimum-norm least-squares
    solution. Returns the ``n x m`` pseudo-inverse.
    """
    U, s, V = svd(A)
    m = len(A)
    n = len(A[0])
    smax = s[0] if s else 0.0
    sinv = [1.0 / sj if sj > rcond * smax else 0.0 for sj in s]
    # A^+ = V diag(sinv) U'
    out = [[0.0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            out[i][j] = sum(V[i][k] * sinv[k] * U[j][k] for k in range(n))
    return out


def condition_number(A):
    """Spectral condition number ``sigma_max / sigma_min`` of ``A``.

    The ratio of the largest to smallest singular value; large means ``A`` is
    ill-conditioned (small perturbations blow up the solution). One for an orthogonal
    matrix; ``inf`` when ``A`` is singular (a zero singular value).
    """
    _, s, _ = svd(A)
    if not s or s[-1] <= 0.0:
        return float("inf")
    return s[0] / s[-1]


def matrix_rank(A, rcond=1e-12):
    """Numerical rank: the number of singular values above ``rcond * sigma_max``.

    Counts the singular directions that carry real signal; singular values below the
    relative tolerance are treated as numerical zeros.
    """
    _, s, _ = svd(A)
    if not s:
        return 0
    smax = s[0]
    if smax <= 0.0:
        return 0
    return sum(1 for sj in s if sj > rcond * smax)


def spectral_norm(A):
    """Spectral (operator 2-) norm: the largest singular value of ``A``."""
    _, s, _ = svd(A)
    return s[0] if s else 0.0


def frobenius_norm(A):
    """Frobenius norm ``sqrt(sum a_ij^2)`` -- equivalently ``sqrt(sum sigma_k^2)``."""
    return sum(a * a for row in A for a in row) ** 0.5
