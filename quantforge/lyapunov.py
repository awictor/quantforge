"""Lyapunov equation solvers (discrete and continuous) and controllability Gramians.

The Lyapunov equation is the linear-algebra core of stability analysis and the LQR/LQG theory.
Discrete: ``A P A^T - P + Q = 0``; continuous: ``A P + P A^T + Q = 0``. For a stable ``A`` and
``Q >= 0`` the solution ``P`` is the system's *Gramian* -- symmetric positive-(semi)definite --
and its existence certifies stability (Lyapunov's theorem). Solving for the *controllability
Gramian* (``Q = B B^T``) quantifies how strongly each state direction can be driven.

Both are linear in the unknown ``P``, so this vectorizes them into an ``n^2 x n^2`` system via the
Kronecker product and solves with :func:`quantforge.lu.lu_solve`, then symmetrizes. Pure standard
library.
"""

from .lu import lu_solve


def _kron(A, B):
    ra, ca = len(A), len(A[0])
    rb, cb = len(B), len(B[0])
    out = [[0.0] * (ca * cb) for _ in range(ra * rb)]
    for i in range(ra):
        for j in range(ca):
            for k in range(rb):
                for l in range(cb):
                    out[i * rb + k][j * cb + l] = A[i][j] * B[k][l]
    return out


def _eye(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def _vec(M):
    # column-major vectorization (stack columns), matching the Kronecker identity
    n = len(M)
    m = len(M[0])
    return [M[i][j] for j in range(m) for i in range(n)]


def _unvec(v, n):
    return [[v[j * n + i] for j in range(n)] for i in range(n)]


def _solve_vectorized(coeff, Q):
    # solve coeff @ vec(P) = -vec(Q) for P (n x n)
    n = len(Q)
    rhs = [-q for q in _vec(Q)]
    p = lu_solve(coeff, rhs)
    P = _unvec(p, n)
    return [[0.5 * (P[i][j] + P[j][i]) for j in range(n)] for i in range(n)]   # symmetrize


def solve_discrete_lyapunov(A, Q):
    """Solve the discrete Lyapunov equation ``A P A^T - P + Q = 0`` for symmetric ``P``.

    Uses ``vec(A P A^T) = (A (x) A) vec(P)``, so ``(A (x) A - I) vec(P) = -vec(Q)``. Requires ``A``
    stable (spectral radius < 1) for a unique solution. Returns the symmetric ``P``.
    """
    n = len(A)
    AkA = _kron(A, A)
    coeff = [[AkA[i][j] - (1.0 if i == j else 0.0) for j in range(n * n)] for i in range(n * n)]
    return _solve_vectorized(coeff, Q)


def solve_continuous_lyapunov(A, Q):
    """Solve the continuous Lyapunov equation ``A P + P A^T + Q = 0`` for symmetric ``P``.

    Uses ``vec(A P + P A^T) = (I (x) A + A (x) I) vec(P)``. Requires ``A`` Hurwitz (eigenvalue real
    parts < 0) for a unique solution. Returns the symmetric ``P``.
    """
    n = len(A)
    I = _eye(n)
    coeff = _kron(I, A)
    AkI = _kron(A, I)
    coeff = [[coeff[i][j] + AkI[i][j] for j in range(n * n)] for i in range(n * n)]
    return _solve_vectorized(coeff, Q)


def controllability_gramian(A, B, discrete=True):
    """Controllability Gramian: Lyapunov solution with ``Q = B B^T``.

    Discrete: solves ``A W A^T - W + B B^T = 0``; continuous: ``A W + W A^T + B B^T = 0``. The
    Gramian ``W`` is positive-definite iff the system is controllable, and its eigen-directions
    rank how easily each mode is driven.
    """
    m = len(B[0])
    n = len(B)
    BBt = [[sum(B[i][k] * B[j][k] for k in range(m)) for j in range(n)] for i in range(n)]
    return solve_discrete_lyapunov(A, BBt) if discrete else solve_continuous_lyapunov(A, BBt)
