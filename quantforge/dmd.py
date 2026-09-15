"""Dynamic mode decomposition (DMD): data-driven spatiotemporal modes.

DMD (Schmid 2010) extracts the dominant oscillation/growth modes of a dynamical system directly
from snapshots -- no governing equations needed. Given states sampled at equal time steps and
stacked as columns of ``X = [x_0, ..., x_{m-1}]`` and ``X' = [x_1, ..., x_m]``, it finds the
best-fit linear operator ``A`` with ``X' ~ A X`` and returns its eigenvalues and eigenvectors
(the *DMD modes*). Each eigenvalue ``lambda_i`` gives a growth rate ``|lambda_i|`` and frequency
``arg(lambda_i)/dt``; the modes are the corresponding spatial patterns. It is the linear-algebra
core of Koopman analysis and reduced-order modelling of fluids, power grids, and video.

This is the exact-DMD variant (Tu et al. 2014): SVD-reduce ``X``, build the compressed operator,
eigendecompose it, and lift the eigenvectors back. Pure standard library on top of
:mod:`quantforge.svd` and :mod:`quantforge.eigen_general`.
"""

from .svd import svd
from .eigen_general import eigenvalues_general


def _matmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][p] * B[p][j] for p in range(k)) for j in range(m)] for i in range(n)]


def _transpose(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def dmd(snapshots, rank=None, dt=1.0):
    """Dynamic mode decomposition of a time-ordered list of state vectors ``snapshots``.

    ``snapshots[t]`` is the state at step ``t`` (a length-``d`` list). Builds ``X`` (columns
    ``0..m-1``) and ``X'`` (columns ``1..m``), SVD-reduces to ``rank`` (default: full numerical
    rank), and returns a dict with ``eigenvalues`` (the discrete-time DMD eigenvalues, possibly
    complex), ``growth_rates`` (``|lambda|``), ``frequencies`` (``arg(lambda)/dt``), and the
    reduced operator ``atilde``.
    """
    m = len(snapshots) - 1
    if m < 1:
        raise ValueError("need at least two snapshots")
    d = len(snapshots[0])
    # X = columns x_0..x_{m-1}  (d x m);  Xp = columns x_1..x_m
    X = [[snapshots[t][i] for t in range(m)] for i in range(d)]
    Xp = [[snapshots[t + 1][i] for t in range(m)] for i in range(d)]

    # SVD of X. quantforge.svd needs rows >= cols; d >= m is the usual tall case.
    if d >= m:
        U, S, V = svd(X)          # U: d x m', S: m', V: m x m'
        Ur, Sr, Vr = U, S, V
    else:
        # wide X: SVD the transpose and swap roles
        Ut, S, Vt = svd(_transpose(X))   # X^T = Ut S Vt^T -> X = Vt S Ut^T
        Ur = Vt                           # d x r
        Sr = S
        Vr = Ut                           # m x r

    r = len([s for s in Sr if s > 1e-12])
    if rank is not None:
        r = min(r, rank)
    Ur = [[Ur[i][j] for j in range(r)] for i in range(d)]
    Vr = [[Vr[i][j] for j in range(r)] for i in range(m)]
    Sinv = [1.0 / Sr[j] for j in range(r)]

    # Atilde = U^T X' V S^-1   (r x r)
    UtXp = _matmul(_transpose(Ur), Xp)          # r x m
    UtXpV = _matmul(UtXp, Vr)                   # r x r
    Atilde = [[UtXpV[i][j] * Sinv[j] for j in range(r)] for i in range(r)]

    eig = eigenvalues_general(Atilde)
    growth = [abs(complex(z)) for z in eig]
    import cmath
    freqs = [cmath.phase(complex(z)) / dt for z in eig]
    return {"eigenvalues": eig, "growth_rates": growth, "frequencies": freqs,
            "atilde": Atilde, "rank": r}
