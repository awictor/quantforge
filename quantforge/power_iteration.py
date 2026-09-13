"""Power iteration and its variants for individual eigenpairs.

Full eigen-decomposition (``jacobi_eigen``) computes every eigenpair; when only one is
needed these are cheaper and only require matrix-vector products:

  * ``power_iteration`` -- converges to the eigenvector of *largest* magnitude
    eigenvalue by repeatedly multiplying by ``A`` and renormalizing; the Rayleigh
    quotient gives the eigenvalue.
  * ``inverse_iteration`` -- power iteration on ``(A - mu I)^{-1}`` converges to the
    eigenvalue *closest* to a shift ``mu`` (e.g. ``mu = 0`` for the smallest-magnitude
    eigenvalue), the standard way to refine a known approximate eigenvalue.
  * ``rayleigh_quotient`` -- ``x' A x / x' x``, the best eigenvalue estimate for a
    given vector.

Pure standard library.
"""

import math

from .portopt import _invert


def _matvec(A, x):
    return [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]


def rayleigh_quotient(A, x):
    """Rayleigh quotient ``x' A x / x' x`` -- the least-squares eigenvalue for ``x``."""
    Ax = _matvec(A, x)
    num = sum(x[i] * Ax[i] for i in range(len(x)))
    den = sum(xi * xi for xi in x)
    if den == 0.0:
        raise ValueError("x must be non-zero")
    return num / den


def _normalize(x):
    nrm = math.sqrt(sum(xi * xi for xi in x))
    if nrm == 0.0:
        raise ValueError("zero vector")
    return [xi / nrm for xi in x]


def power_iteration(A, x0=None, tol=1e-12, max_iter=1000):
    """Dominant eigenpair of ``A`` by power iteration.

    Returns a dict with ``eigenvalue`` (largest magnitude, via the Rayleigh quotient),
    ``eigenvector`` (unit norm), ``n_iter`` and ``converged``. Converges when the
    dominant eigenvalue is unique in magnitude; the sign convention makes the first
    non-negligible component positive.
    """
    n = len(A)
    x = _normalize([1.0] * n if x0 is None else list(x0))
    lam = rayleigh_quotient(A, x)
    n_iter = 0
    converged = False
    for it in range(1, max_iter + 1):
        n_iter = it
        y = _matvec(A, x)
        x_new = _normalize(y)
        lam_new = rayleigh_quotient(A, x_new)
        if abs(lam_new - lam) < tol:
            x, lam = x_new, lam_new
            converged = True
            break
        x, lam = x_new, lam_new
    # Sign convention: first sizeable component positive.
    for xi in x:
        if abs(xi) > 1e-9:
            if xi < 0:
                x = [-v for v in x]
            break
    return {"eigenvalue": lam, "eigenvector": x, "n_iter": n_iter,
            "converged": converged}


def inverse_iteration(A, mu=0.0, x0=None, tol=1e-12, max_iter=1000):
    """Eigenpair of ``A`` whose eigenvalue is closest to the shift ``mu``.

    Runs power iteration on ``(A - mu I)^{-1}``. With ``mu = 0`` this finds the
    smallest-magnitude eigenvalue; with ``mu`` near a known approximate eigenvalue it
    refines that one. Returns the same dict shape as :func:`power_iteration`.
    """
    n = len(A)
    shifted = [[A[i][j] - (mu if i == j else 0.0) for j in range(n)] for i in range(n)]
    inv = _invert(shifted)
    x = _normalize([1.0] * n if x0 is None else list(x0))
    lam = rayleigh_quotient(A, x)
    n_iter = 0
    converged = False
    for it in range(1, max_iter + 1):
        n_iter = it
        y = _matvec(inv, x)
        x_new = _normalize(y)
        lam_new = rayleigh_quotient(A, x_new)
        if abs(lam_new - lam) < tol:
            x, lam = x_new, lam_new
            converged = True
            break
        x, lam = x_new, lam_new
    for xi in x:
        if abs(xi) > 1e-9:
            if xi < 0:
                x = [-v for v in x]
            break
    return {"eigenvalue": lam, "eigenvector": x, "n_iter": n_iter,
            "converged": converged}
