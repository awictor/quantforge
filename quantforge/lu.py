"""LU decomposition with partial pivoting, linear solve, and determinant.

Factors a square matrix as ``P A = L U`` with ``P`` a row-permutation, ``L`` unit
lower-triangular and ``U`` upper-triangular, via Gaussian elimination with partial
pivoting (choosing the largest-magnitude pivot each column for stability). The
factorization solves ``A x = b`` by forward/back substitution and gives the
determinant as the signed product of the ``U`` diagonal. Pure standard library.
"""


def lu_decomposition(A):
    """LU decomposition with partial pivoting of a square matrix.

    Returns ``(L, U, piv, sign)`` where ``L`` is unit lower-triangular, ``U`` upper-
    triangular, ``piv`` the row-permutation (as a list mapping output row -> source
    row) so that ``A[piv] = L U``, and ``sign`` the permutation parity (+1/-1) used
    by :func:`determinant`. Raises on a singular matrix.
    """
    n = len(A)
    if n == 0 or any(len(row) != n for row in A):
        raise ValueError("A must be a non-empty square matrix")
    U = [row[:] for row in A]
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    piv = list(range(n))
    sign = 1

    for k in range(n):
        # Partial pivot: largest |U[i][k]| for i >= k.
        p = max(range(k, n), key=lambda i: abs(U[i][k]))
        if abs(U[p][k]) < 1e-300:
            raise ValueError("matrix is singular")
        if p != k:
            U[k], U[p] = U[p], U[k]
            piv[k], piv[p] = piv[p], piv[k]
            sign = -sign
            # Swap the already-computed part of L.
            for j in range(k):
                L[k][j], L[p][j] = L[p][j], L[k][j]
        for i in range(k + 1, n):
            f = U[i][k] / U[k][k]
            L[i][k] = f
            for j in range(k, n):
                U[i][j] -= f * U[k][j]
    return L, U, piv, sign


def lu_solve(A, b):
    """Solve ``A x = b`` via LU with partial pivoting.

    Factorizes ``A``, permutes ``b``, then forward- and back-substitutes. Returns the
    solution vector. Raises on a singular matrix.
    """
    n = len(A)
    if len(b) != n:
        raise ValueError("b length must match A")
    L, U, piv, _ = lu_decomposition(A)
    pb = [b[piv[i]] for i in range(n)]
    # Forward substitution: L y = pb.
    y = [0.0] * n
    for i in range(n):
        y[i] = pb[i] - sum(L[i][j] * y[j] for j in range(i))
    # Back substitution: U x = y.
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]
    return x


def determinant(A):
    """Determinant of a square matrix as the signed product of the LU pivots.

    ``det = sign * prod(U[i][i])``. Returns 0 for a singular matrix.
    """
    try:
        _, U, _, sign = lu_decomposition(A)
    except ValueError:
        return 0.0
    d = float(sign)
    for i in range(len(U)):
        d *= U[i][i]
    return d
