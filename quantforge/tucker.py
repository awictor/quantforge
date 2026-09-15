"""Tucker decomposition / HOSVD of a 3-way tensor.

Where CP (:mod:`quantforge.cp_decomp`) writes a tensor as a sum of rank-one terms, the Tucker
decomposition factors it into a small *core tensor* multiplied by an orthonormal factor matrix
along each mode:

    X ~ G x_1 U x_2 V x_3 W ,

with ``U (I x r1)``, ``V (J x r2)``, ``W (K x r3)`` having orthonormal columns and core
``G (r1 x r2 x r3)``. It is the tensor analogue of the SVD's ``U S V^T`` and the basis for
multilinear PCA and tensor compression. The higher-order SVD (HOSVD, De Lathauwer et al.) gives
the standard factors: ``U``/``V``/``W`` are the leading left singular vectors of the mode-1/2/3
unfoldings, and the core is the tensor projected onto them.

Full (untruncated) ranks reconstruct ``X`` exactly; truncating the ranks gives a compressed
approximation. Pure standard library on top of :mod:`quantforge.svd`.
"""

from .svd import svd


def _shape(X):
    return len(X), len(X[0]), len(X[0][0])


def _unfold(X, mode):
    I, J, K = _shape(X)
    if mode == 0:
        return [[X[i][j][k] for j in range(J) for k in range(K)] for i in range(I)]
    if mode == 1:
        return [[X[i][j][k] for i in range(I) for k in range(K)] for j in range(J)]
    return [[X[i][j][k] for i in range(I) for j in range(J)] for k in range(K)]


def _transpose(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def _leading_left_vectors(M, r):
    # top-r left singular vectors of M (m x n). svd needs rows >= cols.
    m, n = len(M), len(M[0])
    if m >= n:
        U, S, V = svd(M)
        L = U                       # m x r'
    else:
        # M = (M^T)^T; left vectors of M are right vectors of M^T -> columns of V(M^T)
        Ut, S, Vt = svd(_transpose(M))   # M^T = Ut S Vt^T
        L = Vt                       # m x r'  (right vectors of M^T = left of M)
    r = min(r, len(L[0]))
    return [[L[i][t] for t in range(r)] for i in range(m)]


def _mode_mult(X, Mt, mode):
    # multiply tensor X by matrix Mt (transpose of factor) along `mode`, contracting that index
    I, J, K = _shape(X)
    if mode == 0:
        p = len(Mt)                 # r1
        return [[[sum(Mt[a][i] * X[i][j][k] for i in range(I)) for k in range(K)]
                 for j in range(J)] for a in range(p)]
    if mode == 1:
        p = len(Mt)
        return [[[sum(Mt[b][j] * X[i][j][k] for j in range(J)) for k in range(K)]
                 for b in range(p)] for i in range(len(X))]
    p = len(Mt)
    return [[[sum(Mt[c][k] * X[i][j][k] for k in range(K)) for c in range(p)]
             for j in range(J)] for i in range(len(X))]


def tucker_hosvd(X, ranks=None):
    """Higher-order SVD (Tucker) of a 3-way tensor ``X``.

    ``ranks`` is a triple ``(r1, r2, r3)`` of truncation ranks per mode (default: full ranks).
    Returns a dict with factor matrices ``U``, ``V``, ``W`` (orthonormal columns) and the core
    tensor ``core`` such that ``X ~ core x_1 U x_2 V x_3 W``.
    """
    I, J, K = _shape(X)
    if ranks is None:
        ranks = (I, J, K)
    r1, r2, r3 = ranks
    U = _leading_left_vectors(_unfold(X, 0), r1)
    V = _leading_left_vectors(_unfold(X, 1), r2)
    W = _leading_left_vectors(_unfold(X, 2), r3)
    # core = X x_1 U^T x_2 V^T x_3 W^T
    core = _mode_mult(X, _transpose(U), 0)
    core = _mode_mult(core, _transpose(V), 1)
    core = _mode_mult(core, _transpose(W), 2)
    return {"U": U, "V": V, "W": W, "core": core}


def tucker_reconstruct(U, V, W, core):
    """Reconstruct ``X = core x_1 U x_2 V x_3 W`` from Tucker factors and core."""
    r1 = len(core)
    r2 = len(core[0])
    r3 = len(core[0][0])
    I, J, K = len(U), len(V), len(W)
    out = [[[0.0] * K for _ in range(J)] for _ in range(I)]
    for a in range(r1):
        for b in range(r2):
            for c in range(r3):
                g = core[a][b][c]
                if g == 0.0:
                    continue
                for i in range(I):
                    uia = U[i][a]
                    if uia == 0.0:
                        continue
                    for j in range(J):
                        uv = uia * V[j][b]
                        if uv == 0.0:
                            continue
                        for k in range(K):
                            out[i][j][k] += g * uv * W[k][c]
    return out
