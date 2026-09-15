"""CUR decomposition: interpretable low-rank factorization from actual rows and columns.

SVD/randomized-SVD give the optimal low-rank factors, but their singular vectors are dense linear
combinations of every feature -- hard to interpret. CUR (Mahoney & Drineas 2009) instead picks
``c`` *actual columns* ``C`` and ``r`` *actual rows* ``R`` of the data and a small linking matrix
``U`` so that ``A ~ C U R``. Because ``C`` and ``R`` are real columns/rows, they stay meaningful
(actual genes, documents, assets), which is why CUR is used where interpretability matters.

Columns/rows are chosen by *leverage scores* -- the squared norms of the rows of the top-``k``
right/left singular vectors, which measure how much each column/row contributes to the dominant
subspace. The linking matrix is ``U = C^+ A R^+`` (pseudo-inverses). A seeded
:class:`quantforge.pcg.PCG32` drives the randomized selection. Pure standard library on top of
:mod:`quantforge.svd`.
"""

from .pcg import PCG32
from .svd import svd, pseudo_inverse


def _transpose(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def _matmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][p] * B[p][j] for p in range(k)) for j in range(m)] for i in range(n)]


def column_leverage_scores(A, k):
    """Column leverage scores of ``A`` from its top-``k`` right singular vectors.

    ``score_j = (1/k) sum_{t<k} V[j][t]^2`` -- the normalized contribution of column ``j`` to the
    rank-``k`` dominant subspace. Scores sum to 1 and highlight the most influential columns.
    """
    m, n = len(A), len(A[0])
    # svd needs rows >= cols; SVD A^T if wide so right vectors of A come out as left of A^T
    if m >= n:
        _, S, V = svd(A)                # V: n x r, right singular vectors in columns
        rv = V
    else:
        U, S, _ = svd(_transpose(A))    # A^T = U S W^T -> right vectors of A are columns of U
        rv = U                          # n x r
    r = min(k, len(S))
    return [sum(rv[j][t] ** 2 for t in range(r)) / r for j in range(n)]


def cur(A, c, r, k=None, seed=12345):
    """CUR decomposition ``A ~ C U R`` selecting ``c`` columns and ``r`` rows by leverage scores.

    ``k`` (default ``min(c, r)``) sets the rank used for the leverage scores. Columns and rows are
    sampled without replacement with probability proportional to their leverage. Returns a dict
    with ``C`` (m x c), ``U`` (c x r), ``R`` (r x n), the selected ``col_indices``/``row_indices``,
    and the reconstruction ``error`` (``||A - C U R||_F``).
    """
    m, n = len(A), len(A[0])
    if not (1 <= c <= n and 1 <= r <= m):
        raise ValueError("require 1 <= c <= ncols and 1 <= r <= nrows")
    if k is None:
        k = min(c, r)
    rng = PCG32(seed)

    col_lev = column_leverage_scores(A, k)
    row_lev = column_leverage_scores(_transpose(A), k)

    col_idx = sorted(_sample_without_replacement(col_lev, c, rng))
    row_idx = sorted(_sample_without_replacement(row_lev, r, rng))

    C = [[A[i][j] for j in col_idx] for i in range(m)]
    R = [[A[i][j] for j in range(n)] for i in row_idx]
    # U = C^+ A R^+. quantforge.pseudo_inverse (via svd) needs rows >= cols, so for a wide
    # matrix M compute M^+ = (M^T ^+)^T which keeps the SVD tall-skinny.
    def _pinv(M):
        if len(M) >= len(M[0]):
            return pseudo_inverse(M)
        return _transpose(pseudo_inverse(_transpose(M)))

    Cp = _pinv(C)                      # c x m
    Rp = _pinv(R)                      # n x r
    U = _matmul(_matmul(Cp, A), Rp)    # c x r

    CUR = _matmul(_matmul(C, U), R)
    err = sum((A[i][j] - CUR[i][j]) ** 2 for i in range(m) for j in range(n)) ** 0.5
    return {"C": C, "U": U, "R": R, "col_indices": col_idx, "row_indices": row_idx,
            "error": err}


def _sample_without_replacement(scores, count, rng):
    # probability-proportional-to-score sampling without replacement
    idx = list(range(len(scores)))
    w = list(scores)
    chosen = []
    for _ in range(min(count, len(idx))):
        total = sum(w[i] for i in idx)
        if total <= 0:
            chosen.extend(idx[:count - len(chosen)])
            break
        target = rng.random() * total
        acc = 0.0
        pick = idx[-1]
        for i in idx:
            acc += w[i]
            if acc >= target:
                pick = i
                break
        chosen.append(pick)
        idx.remove(pick)
    return chosen
