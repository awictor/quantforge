"""Matrix-pencil method: SVD-based estimation of exponential modes.

Prony's method (:mod:`quantforge.prony`) recovers a sum of exponentials from the *minimum* number
of samples, but is fragile under noise because it roots a polynomial fit to a single recurrence.
The matrix-pencil method (Hua & Sarkar 1990) is the noise-robust alternative: it stacks the data
into a Hankel matrix, truncates it to rank ``p`` with an SVD (discarding the noise subspace), and
recovers the modes ``z_i`` as the eigenvalues of a pencil built from the truncated matrix. This
is the ESPRIT-style estimator behind modal analysis and direction-of-arrival.

Given samples ``y_k`` modelled as ``sum_i a_i z_i^k``, returns the modes ``z_i`` (and, via
:func:`quantforge.prony`-style least squares, the amplitudes). Pure standard library on top of
:mod:`quantforge.svd` and :mod:`quantforge.eigen_general`.
"""

from .svd import svd
from .eigen_general import eigenvalues_general
from .prony import _lstsq_complex


def _matmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][p] * B[p][j] for p in range(k)) for j in range(m)] for i in range(n)]


def _transpose(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def _pinv_from_svd(A, rank):
    # Moore-Penrose pseudo-inverse truncated to `rank` singular values
    U, S, V = svd(A)
    r = min(rank, len(S))
    # A+ = V (1/S) U^T, using only the top r components
    m, n = len(A), len(A[0])
    Ap = [[0.0] * m for _ in range(n)]
    for t in range(r):
        if S[t] < 1e-14:
            continue
        inv_s = 1.0 / S[t]
        for i in range(n):
            vi = V[i][t]
            for j in range(m):
                Ap[i][j] += vi * inv_s * U[j][t]
    return Ap


def matrix_pencil(y, p, pencil=None):
    """Estimate ``p`` exponential modes of the sample sequence ``y`` by the matrix-pencil method.

    ``pencil`` is the pencil parameter ``L`` (Hankel column count minus one); if ``None`` it
    defaults to ``len(y) // 3``, near the noise-optimal choice. Requires ``p <= L <= len(y) - p``.
    Returns ``(modes, amplitudes)`` with ``y_k ~ sum amplitudes[i] * modes[i]**k``; both complex.
    """
    N = len(y)
    if p < 1:
        raise ValueError("need p >= 1")
    L = pencil if pencil is not None else N // 3
    if not (p <= L <= N - p):
        raise ValueError("require p <= L <= N - p")
    yc = [complex(v) for v in y]

    # Hankel data matrix Y0 with rows k=0..N-L-1, cols j=0..L
    rows = N - L
    Y = [[yc[k + j] for j in range(L + 1)] for k in range(rows)]
    # Y1 = Y[:, :L], Y2 = Y[:, 1:]  (shift by one column)
    Y1 = [[Y[i][j] for j in range(L)] for i in range(rows)]
    Y2 = [[Y[i][j + 1] for j in range(L)] for i in range(rows)]

    # SVD works on real matrices here; split into real/imag if data is complex.
    # For real input (the common case) Y1/Y2 are real, so use the real SVD directly.
    Y1r = [[v.real for v in row] for row in Y1]
    Y1p = _pinv_from_svd(Y1r, p)                # p x rows, real
    # pencil matrix M = Y1^+ Y2 (rows: L x L), truncated to rank p implicitly by the pinv
    Y2r = [[v.real for v in row] for row in Y2]
    M = _matmul(Y1p, Y2r)                       # L x L real
    modes = eigenvalues_general(M)
    # keep the p largest-magnitude eigenvalues (signal subspace), drop near-zero noise modes
    modes = sorted(modes, key=lambda z: -abs(z))[:p]

    # amplitudes by least squares on the Vandermonde system
    V = [[modes[i] ** k for i in range(p)] for k in range(N)]
    amps = _lstsq_complex(V, yc)
    return modes, amps
