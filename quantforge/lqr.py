"""Discrete linear-quadratic regulator (LQR) and controllability.

The LQR is the foundational optimal-control result: for a linear system ``x_{k+1} = A x_k + B u_k``
and a quadratic cost ``sum x^T Q x + u^T R u``, the cost-minimizing control is a *linear state
feedback* ``u_k = -K x_k``. The gain ``K`` comes from the stabilizing solution ``P`` of the
discrete algebraic Riccati equation (DARE)

    P = A^T P A - A^T P B (R + B^T P B)^{-1} B^T P A + Q ,

which this solves by fixed-point iteration. The resulting closed-loop system ``A - B K`` is
stable (all eigenvalues inside the unit circle) whenever the system is controllable. Also
provides the :func:`controllability_matrix` and a rank-based controllability test. Pure standard
library on top of :func:`quantforge.lu.lu_solve`.
"""

from .lu import lu_solve


def _matmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][p] * B[p][j] for p in range(k)) for j in range(m)] for i in range(n)]


def _transpose(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def _add(A, B, sign=1.0):
    return [[A[i][j] + sign * B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def _inv(M):
    n = len(M)
    cols = [lu_solve(M, [1.0 if i == j else 0.0 for i in range(n)]) for j in range(n)]
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def dare(A, B, Q, R, tol=1e-12, max_iter=1000):
    """Stabilizing solution ``P`` of the discrete algebraic Riccati equation, by iteration.

    ``A`` (n x n), ``B`` (n x m), ``Q`` (n x n, >=0), ``R`` (m x m, >0). Returns the symmetric
    ``P`` that the DARE fixed point converges to. Raises if it fails to converge.
    """
    n = len(A)
    At = _transpose(A)
    Bt = _transpose(B)
    P = [[Q[i][j] for j in range(n)] for i in range(n)]
    for _ in range(max_iter):
        # S = R + B^T P B  (m x m)
        BtP = _matmul(Bt, P)
        S = _add(R, _matmul(BtP, B))
        # K_term = A^T P B S^{-1} B^T P A
        AtP = _matmul(At, P)
        AtPB = _matmul(AtP, B)
        Sinv = _inv(S)
        term = _matmul(_matmul(AtPB, Sinv), _matmul(BtP, A))
        P_new = _add(_add(_matmul(AtP, A), term, sign=-1.0), Q)
        # symmetrize to counter drift
        P_new = [[0.5 * (P_new[i][j] + P_new[j][i]) for j in range(n)] for i in range(n)]
        diff = max(abs(P_new[i][j] - P[i][j]) for i in range(n) for j in range(n))
        P = P_new
        if diff < tol:
            return P
    return P


def lqr(A, B, Q, R):
    """Discrete LQR gain ``K`` and Riccati solution ``P`` for ``x_{k+1}=Ax+Bu``, cost ``x'Qx+u'Ru``.

    Returns a dict with ``K`` (m x n feedback gain; optimal control is ``u = -K x``) and ``P``
    (the DARE solution). The closed-loop ``A - B K`` is stable when the pair ``(A, B)`` is
    controllable.
    """
    P = dare(A, B, Q, R)
    Bt = _transpose(B)
    BtP = _matmul(Bt, P)
    S = _add(R, _matmul(BtP, B))                      # R + B^T P B
    K = _matmul(_inv(S), _matmul(BtP, A))             # (R + B^T P B)^{-1} B^T P A
    return {"K": K, "P": P}


def controllability_matrix(A, B):
    """Controllability matrix ``[B, AB, A^2 B, ..., A^{n-1} B]`` (n x n*m)."""
    n = len(A)
    m = len(B[0])
    blocks = [B]
    AkB = B
    for _ in range(1, n):
        AkB = _matmul(A, AkB)
        blocks.append(AkB)
    # horizontally stack
    return [[blocks[b][i][j] for b in range(n) for j in range(m)] for i in range(n)]


def is_controllable(A, B, tol=1e-9):
    """True if ``(A, B)`` is controllable (controllability matrix has full row rank ``n``)."""
    from .svd import matrix_rank
    C = controllability_matrix(A, B)
    return matrix_rank(C, tol) == len(A)
