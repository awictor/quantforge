"""Balanced truncation model-order reduction via Hankel singular values.

A high-order LTI system often has states that barely affect the input-output behaviour --
weakly controllable *and* weakly observable. Balanced truncation (Moore 1981) identifies and
removes them: it balances the realization so the controllability and observability Gramians are
equal and diagonal, then discards the states with the smallest *Hankel singular values* (the
shared diagonal). Those HSVs measure each state's input-output energy, are invariant under
coordinate change, and give a hard error bound on the reduced model. It is the standard
principled way to shrink a model while keeping its dynamics.

:func:`hankel_singular_values` are the sqrt-eigenvalues of ``Wc Wo`` (product of the two
Gramians); :func:`balanced_truncation` returns a reduced ``(A, B, C, D)`` keeping the top ``r``.
Pure standard library on top of :mod:`quantforge.lyapunov` and :func:`quantforge.pca.jacobi_eigen`.
"""

import math

from .lyapunov import solve_discrete_lyapunov
from .pca import jacobi_eigen


def _matmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][p] * B[p][j] for p in range(k)) for j in range(m)] for i in range(n)]


def _transpose(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def controllability_gramian_d(A, B):
    """Discrete controllability Gramian ``Wc``: solves ``A Wc A^T - Wc + B B^T = 0``."""
    n = len(A)
    m = len(B[0])
    BBt = [[sum(B[i][k] * B[j][k] for k in range(m)) for j in range(n)] for i in range(n)]
    return solve_discrete_lyapunov(A, BBt)


def observability_gramian_d(A, C):
    """Discrete observability Gramian ``Wo``: solves ``A^T Wo A - Wo + C^T C = 0``."""
    n = len(A)
    p = len(C)
    CtC = [[sum(C[k][i] * C[k][j] for k in range(p)) for j in range(n)] for i in range(n)]
    return solve_discrete_lyapunov(_transpose(A), CtC)


def hankel_singular_values(A, B, C):
    """Hankel singular values of a stable discrete system: ``sqrt(eig(Wc Wo))``, descending.

    They quantify each state's input-output energy and are invariant under state-coordinate
    change; the smallest ones flag states that can be truncated with little effect.
    """
    Wc = controllability_gramian_d(A, B)
    Wo = observability_gramian_d(A, C)
    WcWo = _matmul(Wc, Wo)
    # WcWo is not symmetric, but its eigenvalues are real and non-negative (product of two SPD).
    # Symmetrize the similar matrix via Wo-weighting is complex; here use eig of the symmetric
    # part is wrong. Instead form the SPD product's eigenvalues through Wc^{1/2} Wo Wc^{1/2}.
    # Simpler and exact for the HSVs: eigenvalues of Wc Wo equal eigenvalues of the symmetric
    # M = L^T Wo L where Wc = L L^T. Use a Cholesky-free route via jacobi on (Wc Wo + (Wc Wo)^T)/2
    # only when it is already near-symmetric; for robustness compute L by eigen-sqrt of Wc.
    vals, vecs = jacobi_eigen(Wc)                  # Wc = V diag(vals) V^T, vecs are ROWS
    n = len(A)
    # Wc^{1/2} = V diag(sqrt vals) V^T
    sq = [max(v, 0.0) ** 0.5 for v in vals]
    Whalf = [[sum(vecs[t][i] * sq[t] * vecs[t][j] for t in range(n)) for j in range(n)]
             for i in range(n)]
    M = _matmul(_matmul(Whalf, Wo), Whalf)         # symmetric SPD, similar to Wc Wo
    ev, _ = jacobi_eigen(M)
    hsv = sorted((max(e, 0.0) ** 0.5 for e in ev), reverse=True)
    return hsv


def balanced_truncation(A, B, C, D, r):
    """Reduce a stable discrete LTI ``(A, B, C, D)`` to order ``r`` by balanced truncation.

    Balances the realization (so ``Wc = Wo = diag(HSV)``) and keeps the ``r`` states with the
    largest Hankel singular values. Returns the reduced ``(Ar, Br, Cr, Dr)``. ``D`` is unchanged.
    """
    n = len(A)
    if not (1 <= r <= n):
        raise ValueError("require 1 <= r <= n")
    Wc = controllability_gramian_d(A, B)
    Wo = observability_gramian_d(A, C)
    # balancing transform T from the Gramians (square-root method)
    vc, Vc = jacobi_eigen(Wc)                      # rows are eigenvectors
    sq = [max(v, 0.0) ** 0.5 for v in vc]
    # Lc = Vc^T diag(sqrt) : columns; build L such that Wc = L L^T
    L = [[Vc[t][i] * sq[t] for t in range(n)] for i in range(n)]     # n x n, L[i][t]
    # M = L^T Wo L
    Lt = _transpose(L)
    M = _matmul(_matmul(Lt, Wo), L)
    sv, U = jacobi_eigen(M)                          # M = U diag(sv) U^T (rows)
    order = sorted(range(n), key=lambda t: -sv[t])
    hsv = [max(sv[t], 0.0) ** 0.5 for t in order]
    # balancing transform T = L U diag(sv^{-1/4}), T^{-1} = diag(sv^{1/4}) U^T L^{-1}
    # build columns of T in sorted order
    Ucols = [[U[order[t]][i] for i in range(n)] for t in range(n)]   # Ucols[t] = t-th eigvec
    sinv4 = [max(sv[order[t]], 1e-30) ** -0.25 for t in range(n)]
    s4 = [max(sv[order[t]], 1e-30) ** 0.25 for t in range(n)]
    # T[:,t] = (L U[:,t]) * sinv4[t]
    LU = [[sum(L[i][k] * Ucols[t][k] for k in range(n)) for t in range(n)] for i in range(n)]
    T = [[LU[i][t] * sinv4[t] for t in range(n)] for i in range(n)]
    # T^{-1} row t = s4[t] * (U[:,t]^T L^{-1}); compute via solving, but easier: Tinv = diag(s4) U^T L^{-1}
    from .lu import lu_solve
    Linv_cols = [lu_solve(L, [1.0 if i == j else 0.0 for i in range(n)]) for j in range(n)]
    Linv = [[Linv_cols[j][i] for j in range(n)] for i in range(n)]
    UtLinv = [[sum(Ucols[t][k] * Linv[k][j] for k in range(n)) for j in range(n)]
              for t in range(n)]
    Tinv = [[s4[t] * UtLinv[t][j] for j in range(n)] for t in range(n)]
    # balanced (Ab, Bb, Cb) = (Tinv A T, Tinv B, C T)
    Ab = _matmul(_matmul(Tinv, A), T)
    Bb = _matmul(Tinv, B)
    Cb = _matmul(C, T)
    # truncate to top r states
    Ar = [[Ab[i][j] for j in range(r)] for i in range(r)]
    Br = [[Bb[i][j] for j in range(len(B[0]))] for i in range(r)]
    Cr = [[Cb[i][j] for j in range(r)] for i in range(len(C))]
    return Ar, Br, Cr, [row[:] for row in D]
