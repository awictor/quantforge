"""CP (CANDECOMP/PARAFAC) decomposition of a 3-way tensor by alternating least squares.

The SVD factors a matrix; the CP decomposition is its tensor analogue. A 3-way tensor
``X[i,j,k]`` of shape ``(I, J, K)`` is approximated by a sum of ``R`` rank-one terms,

    X[i,j,k] ~ sum_{r=1}^{R} A[i,r] B[j,r] C[k,r] ,

with factor matrices ``A (I x R)``, ``B (J x R)``, ``C (K x R)``. Unlike the matrix SVD the CP
decomposition is (generically) *unique* up to scaling and permutation, which is why it recovers
interpretable latent factors in chemometrics, neuroimaging, and multi-way data mining.

This fits it by alternating least squares: hold two factor matrices fixed and solve the third as
a linear least-squares problem (via the Khatri-Rao product and the normal equations), cycling
until the fit stops improving. A seeded :class:`quantforge.pcg.PCG32` initializes the factors.
Pure standard library.
"""

from .pcg import PCG32
from .lu import lu_solve


def _unfold(X, mode):
    # matricize a 3-way tensor X (list of I of J-by-K) along the given mode
    I = len(X)
    J = len(X[0])
    K = len(X[0][0])
    if mode == 0:
        return [[X[i][j][k] for j in range(J) for k in range(K)] for i in range(I)]
    if mode == 1:
        return [[X[i][j][k] for i in range(I) for k in range(K)] for j in range(J)]
    return [[X[i][j][k] for i in range(I) for j in range(J)] for k in range(K)]


def _khatri_rao(P, Q):
    # column-wise Khatri-Rao (columns must match): rows = rows(P)*rows(Q)
    rp, rq = len(P), len(Q)
    R = len(P[0])
    return [[P[a][r] * Q[b][r] for r in range(R)] for a in range(rp) for b in range(rq)]


def _matmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][p] * B[p][j] for p in range(k)) for j in range(m)] for i in range(n)]


def _transpose(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def _solve_factor(unfolded, kr):
    # solve  factor * kr^T ~ unfolded  =>  factor = unfolded * kr * (kr^T kr)^-1
    # normal equations on the R x R Gram matrix (small, well-posed)
    R = len(kr[0])
    G = _matmul(_transpose(kr), kr)             # R x R
    rhs = _matmul(unfolded, kr)                 # (dim) x R
    # solve G^T f^T = rhs^T column by column  ->  factor row f solves f G = rhs_row
    Gt = _transpose(G)
    out = []
    for row in rhs:
        out.append(lu_solve(Gt, row))
    return out                                  # (dim) x R


def cp_als(X, rank, iterations=200, seed=12345, tol=1e-10):
    """CP/PARAFAC decomposition of a 3-way tensor ``X`` into ``rank`` rank-one terms by ALS.

    ``X`` is a nested list of shape ``(I, J, K)``. Returns a dict with factor matrices ``A``, ``B``,
    ``C`` (each ``dim x rank``), the reconstruction ``error`` (Frobenius) and ``n_iter``.
    """
    I, J, K = len(X), len(X[0]), len(X[0][0])
    if rank < 1:
        raise ValueError("rank must be >= 1")
    rng = PCG32(seed)
    B = [[rng.random() for _ in range(rank)] for _ in range(J)]
    C = [[rng.random() for _ in range(rank)] for _ in range(K)]
    A = [[rng.random() for _ in range(rank)] for _ in range(I)]

    X0, X1, X2 = _unfold(X, 0), _unfold(X, 1), _unfold(X, 2)

    def recon_err():
        e = 0.0
        for i in range(I):
            for j in range(J):
                for k in range(K):
                    v = sum(A[i][r] * B[j][r] * C[k][r] for r in range(rank))
                    e += (X[i][j][k] - v) ** 2
        return e ** 0.5

    prev = recon_err()
    n_iter = 0
    for it in range(iterations):
        n_iter = it + 1
        A = _solve_factor(X0, _khatri_rao(B, C))    # I x R, KR over (J,K)
        B = _solve_factor(X1, _khatri_rao(A, C))    # J x R, KR over (I,K)
        C = _solve_factor(X2, _khatri_rao(A, B))    # K x R, KR over (I,J)
        err = recon_err()
        if abs(prev - err) < tol * max(1.0, prev):
            prev = err
            break
        prev = err
    return {"A": A, "B": B, "C": C, "error": prev, "n_iter": n_iter}


def cp_reconstruct(A, B, C):
    """Reconstruct the tensor ``X[i,j,k] = sum_r A[i,r] B[j,r] C[k,r]`` from CP factors."""
    I, J, K, R = len(A), len(B), len(C), len(A[0])
    return [[[sum(A[i][r] * B[j][r] * C[k][r] for r in range(R)) for k in range(K)]
             for j in range(J)] for i in range(I)]
