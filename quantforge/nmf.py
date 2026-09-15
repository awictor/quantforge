"""Non-negative matrix factorization (Lee-Seung multiplicative updates).

NMF factors a non-negative matrix ``V (m x n)`` into two non-negative factors ``W (m x k)`` and
``H (k x n)`` with ``V ~ W H``. The non-negativity makes the factors *additive* and
interpretable -- parts-based representations of images, topics in a document-term matrix,
spectra in a mixture -- where the SVD's signed components are not. This uses the Lee-Seung (2001)
multiplicative update rules, which minimize the Frobenius reconstruction error while keeping both
factors non-negative at every step (no projection needed).

A seeded :class:`quantforge.pcg.PCG32` initializes the factors so runs are reproducible. Pure
standard library.
"""

from .pcg import PCG32


def _matmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][p] * B[p][j] for p in range(k)) for j in range(m)] for i in range(n)]


def _transpose(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def nmf(V, k, iterations=300, seed=12345, tol=1e-9, eps=1e-12):
    """Factor a non-negative matrix ``V`` as ``W H`` with ``W, H >= 0`` and inner dimension ``k``.

    ``V`` is a list of non-negative rows. Runs Lee-Seung multiplicative updates for up to
    ``iterations`` (stopping early if the Frobenius error stalls within ``tol``). Returns a dict
    with ``W`` (m x k), ``H`` (k x n), ``error`` (``||V - WH||_F``) and ``n_iter``.
    """
    m = len(V)
    n = len(V[0])
    if any(V[i][j] < 0 for i in range(m) for j in range(n)):
        raise ValueError("V must be non-negative")
    if not (1 <= k <= min(m, n) + max(m, n)):
        raise ValueError("k out of range")
    rng = PCG32(seed)
    # non-negative random init scaled to the data magnitude
    scale = (sum(V[i][j] for i in range(m) for j in range(n)) / (m * n)) ** 0.5 or 1.0
    W = [[rng.random() * scale for _ in range(k)] for _ in range(m)]
    H = [[rng.random() * scale for _ in range(n)] for _ in range(k)]

    def fro_err(W, H):
        WH = _matmul(W, H)
        return sum((V[i][j] - WH[i][j]) ** 2 for i in range(m) for j in range(n)) ** 0.5

    prev = fro_err(W, H)
    n_iter = 0
    for it in range(iterations):
        n_iter = it + 1
        # H <- H * (W^T V) / (W^T W H)
        Wt = _transpose(W)
        WtV = _matmul(Wt, V)
        WtWH = _matmul(_matmul(Wt, W), H)
        H = [[H[a][j] * WtV[a][j] / (WtWH[a][j] + eps) for j in range(n)] for a in range(k)]
        # W <- W * (V H^T) / (W H H^T)
        Ht = _transpose(H)
        VHt = _matmul(V, Ht)
        WHHt = _matmul(W, _matmul(H, Ht))
        W = [[W[i][a] * VHt[i][a] / (WHHt[i][a] + eps) for a in range(k)] for i in range(m)]

        err = fro_err(W, H)
        if abs(prev - err) < tol * max(1.0, prev):
            prev = err
            break
        prev = err
    return {"W": W, "H": H, "error": prev, "n_iter": n_iter}
