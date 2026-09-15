"""Pole placement by Ackermann's formula.

Where the LQR (:mod:`quantforge.lqr`) chooses the feedback gain to minimize a cost, *pole
placement* lets you specify the closed-loop dynamics directly: given desired eigenvalues for
``A - B K``, find the gain ``K`` that puts them there. For a controllable single-input system,
Ackermann's formula gives ``K`` in closed form,

    K = [0, ..., 0, 1] C^{-1} phi(A) ,

where ``C = [B, AB, ..., A^{n-1}B]`` is the controllability matrix and ``phi`` is the desired
characteristic polynomial ``prod (s - lambda_i)`` evaluated at the matrix ``A``. Placing the poles
sets the transient response (speed, damping) exactly. Pure standard library on top of
:func:`quantforge.lqr.controllability_matrix` and :func:`quantforge.lu.lu_solve`.
"""

from .lqr import controllability_matrix
from .lu import lu_solve


def _matmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][p] * B[p][j] for p in range(k)) for j in range(m)] for i in range(n)]


def _eye(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def _poly_from_roots(roots):
    # expand prod (s - r) into coefficients, highest degree first (leading 1)
    coeffs = [1.0]
    for r in roots:
        new = [0.0] * (len(coeffs) + 1)
        for i, c in enumerate(coeffs):
            new[i] += c
            new[i + 1] += -r * c
        coeffs = new
    return coeffs


def _matrix_poly(A, coeffs):
    # evaluate polynomial (highest-degree-first coeffs) at matrix A via Horner
    n = len(A)
    result = [[coeffs[0] if i == j else 0.0 for j in range(n)] for i in range(n)]
    for c in coeffs[1:]:
        result = _matmul(result, A)
        for i in range(n):
            result[i][i] += c
    return result


def ackermann(A, B, desired_poles):
    """Single-input pole-placement gain ``K`` via Ackermann's formula.

    ``A`` (n x n), ``B`` (n x 1); ``desired_poles`` a list of ``n`` target closed-loop
    eigenvalues (real here). Returns the ``1 x n`` gain ``K`` such that ``A - B K`` has exactly
    those eigenvalues. Requires ``(A, B)`` controllable.
    """
    n = len(A)
    if len(B[0]) != 1:
        raise ValueError("Ackermann's formula is for single-input systems (B is n x 1)")
    if len(desired_poles) != n:
        raise ValueError("need exactly n desired poles")
    C = controllability_matrix(A, B)              # n x n for single input
    phi = _matrix_poly(A, _poly_from_roots(desired_poles))   # phi(A), n x n
    # last row of C^{-1}: solve C^T y = e_n  =>  y^T = e_n^T C^{-1}
    Ct = [[C[j][i] for j in range(n)] for i in range(n)]
    e_n = [0.0] * n
    e_n[-1] = 1.0
    last_row = lu_solve(Ct, e_n)                  # this is e_n^T C^{-1}
    # K = last_row @ phi(A)  (1 x n)
    K = [[sum(last_row[k] * phi[k][j] for k in range(n)) for j in range(n)]]
    return K
