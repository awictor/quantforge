"""Randomized SVD: fast low-rank approximation of a matrix.

The full SVD (:func:`quantforge.svd.svd`) costs ``O(mn min(m,n))`` -- wasteful when you only want
the top ``k`` singular triplets of a large, approximately low-rank matrix. The randomized SVD
(Halko, Martinsson & Tropp 2011) sketches the range of ``A`` with a random Gaussian projection,
orthonormalizes it, and computes an exact SVD of the small projected matrix. A few *power
iterations* sharpen the sketch when the singular values decay slowly. The result approximates the
top-``k`` SVD to high accuracy at a fraction of the cost -- the workhorse behind large-scale PCA,
image compression, and recommender factorizations.

Returns ``(U, s, V)`` with ``U`` (``m x k``), ``s`` (length ``k``, descending) and ``V``
(``n x k``) such that ``A ~ U diag(s) V^T``. A seeded :class:`quantforge.pcg.PCG32` drives the
random projection. Pure standard library on top of :mod:`quantforge.qr` and :mod:`quantforge.svd`.
"""

from .pcg import PCG32
from .particle_filter import pcg_gaussian
from .qr import qr_decomposition
from .svd import svd


def _matmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][p] * B[p][j] for p in range(k)) for j in range(m)] for i in range(n)]


def _transpose(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def randomized_svd(A, k, n_oversample=5, n_power=2, seed=12345):
    """Rank-``k`` randomized SVD of matrix ``A`` (list of rows).

    ``n_oversample`` extra random columns improve accuracy (Halko et al. recommend ~5-10);
    ``n_power`` power iterations sharpen the range estimate for slowly-decaying spectra. Returns
    ``(U, s, V)`` with the leading ``k`` singular triplets, ``A ~ U diag(s) V^T``.
    """
    m = len(A)
    n = len(A[0])
    r = min(k + n_oversample, n, m)
    rng = PCG32(seed)

    # random Gaussian test matrix Omega (n x r) and sketch Y = A Omega (m x r)
    Omega = [[pcg_gaussian(rng, 0.0, 1.0) for _ in range(r)] for _ in range(n)]
    Y = _matmul(A, Omega)

    At = _transpose(A)
    # power iterations: Y <- A (A^T Y), re-orthonormalizing each step for stability
    for _ in range(n_power):
        Q, _ = qr_decomposition(Y)
        Qthin = [[Q[i][j] for j in range(r)] for i in range(m)]
        Z = _matmul(At, Qthin)                 # n x r
        Qz, _ = qr_decomposition(Z)
        Qzthin = [[Qz[i][j] for j in range(r)] for i in range(n)]
        Y = _matmul(A, Qzthin)                 # m x r

    # orthonormal basis Q for the range of Y
    Q, _ = qr_decomposition(Y)
    Qthin = [[Q[i][j] for j in range(r)] for i in range(m)]   # m x r

    # project: B = Q^T A  (r x n). quantforge.svd needs rows >= cols, so SVD B^T (n x r):
    #   B^T = Uc diag(s) Vc^T  =>  B = Vc diag(s) Uc^T, i.e. left vecs of B are Vc.
    B = _matmul(_transpose(Qthin), A)          # r x n
    Bt = _transpose(B)                          # n x r  (n >= r)
    Uc, s, Vc = svd(Bt)                         # Uc: n x r', s: r', Vc: r x r'
    # left singular vectors of B are Vc (r x r'); right singular vectors are Uc (n x r')
    # U = Q * Vc ; keep top k
    U_full = _matmul(Qthin, Vc)                 # m x r'
    kk = min(k, len(s))
    U = [[U_full[i][j] for j in range(kk)] for i in range(m)]
    s_out = s[:kk]
    V = [[Uc[i][j] for j in range(kk)] for i in range(n)]
    return U, s_out, V
