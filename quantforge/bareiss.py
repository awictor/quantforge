"""Exact integer/rational linear algebra (Bareiss fraction-free elimination).

Floating Gaussian elimination accumulates round-off; for integer or rational matrices the
answer can be computed *exactly*. Bareiss's fraction-free algorithm runs Gaussian
elimination while keeping every intermediate entry an exact integer (each division is
guaranteed exact), yielding the determinant with no rounding at all. Paired with exact
``Fraction`` arithmetic it also solves integer linear systems and inverts integer
matrices exactly. These are the tools for exact rank, Cramer-style solutions, and
verifying a floating solver. Pure standard library (uses :mod:`fractions`).
"""

from fractions import Fraction


def bareiss_determinant(matrix):
    """Exact determinant of a square integer (or rational) matrix by Bareiss elimination.

    Fraction-free Gaussian elimination: each step divides by the previous pivot and the
    division is always exact, so with integer input the result is an exact integer (no
    round-off, no overflow on Python big integers). Returns an ``int`` for integer input.
    """
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("matrix must be square and non-empty")
    a = [[Fraction(v) for v in row] for row in matrix]
    sign = 1
    prev = Fraction(1)
    for k in range(n - 1):
        if a[k][k] == 0:
            # Pivot: swap with a lower row that has a non-zero entry in column k.
            swap = None
            for i in range(k + 1, n):
                if a[i][k] != 0:
                    swap = i
                    break
            if swap is None:
                return 0            # singular
            a[k], a[swap] = a[swap], a[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) / prev
            a[i][k] = Fraction(0)
        prev = a[k][k]
    det = sign * a[n - 1][n - 1]
    return int(det) if det.denominator == 1 else det


def rational_solve(A, b):
    """Exact solution ``x`` of ``A x = b`` over the rationals (Gaussian elimination).

    ``A`` is a square integer/rational matrix and ``b`` a right-hand-side vector. Returns
    the exact solution as a list of :class:`fractions.Fraction`. Raises if ``A`` is
    singular. No round-off: an ill-conditioned system that a float solver botches comes
    out exact.
    """
    n = len(A)
    if n == 0 or any(len(row) != n for row in A):
        raise ValueError("A must be square and non-empty")
    if len(b) != n:
        raise ValueError("b length must match A")
    # Augmented matrix in exact rationals.
    m = [[Fraction(A[i][j]) for j in range(n)] + [Fraction(b[i])] for i in range(n)]
    for k in range(n):
        # Partial pivot on the first non-zero entry (exact, so any non-zero works).
        if m[k][k] == 0:
            swap = next((i for i in range(k + 1, n) if m[i][k] != 0), None)
            if swap is None:
                raise ValueError("matrix is singular")
            m[k], m[swap] = m[swap], m[k]
        pivot = m[k][k]
        for i in range(n):
            if i != k and m[i][k] != 0:
                factor = m[i][k] / pivot
                for j in range(k, n + 1):
                    m[i][j] -= factor * m[k][j]
    return [m[i][n] / m[i][i] for i in range(n)]


def rational_inverse(A):
    """Exact inverse of a square integer/rational matrix as a matrix of ``Fraction``.

    Solves ``A X = I`` column by column with :func:`rational_solve`. Raises if ``A`` is
    singular. Multiplying the result by ``A`` returns the exact identity.
    """
    n = len(A)
    if n == 0 or any(len(row) != n for row in A):
        raise ValueError("A must be square and non-empty")
    cols = []
    for c in range(n):
        e = [1 if i == c else 0 for i in range(n)]
        cols.append(rational_solve(A, e))
    # cols[c] is column c of the inverse; transpose to row-major.
    return [[cols[c][r] for c in range(n)] for r in range(n)]
