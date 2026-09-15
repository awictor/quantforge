"""Preconditioners and preconditioned conjugate gradient.

An iterative solver's convergence speed is governed by the condition number of ``A``: the more
spread-out the spectrum, the more iterations. A *preconditioner* ``M ~ A`` (cheap to invert)
transforms the system so ``M^{-1} A`` is better conditioned, cutting iteration counts
dramatically on stiff problems. Two classic choices:

* :func:`jacobi_preconditioner` -- diagonal (Jacobi) preconditioning ``M = diag(A)``; trivial and
  effective when ``A`` is diagonally dominant.
* :func:`incomplete_cholesky` -- IC(0): a Cholesky factor with the *same sparsity pattern* as
  ``A`` (no fill-in), the workhorse preconditioner for symmetric positive-definite systems.

:func:`preconditioned_cg` runs conjugate gradient with an ``M^{-1}`` apply, converging in far
fewer iterations than plain CG on ill-conditioned SPD systems. Pure standard library.
"""

import math


def _dot(a, b):
    return sum(a[i] * b[i] for i in range(len(a)))


def jacobi_preconditioner(A):
    """Diagonal (Jacobi) preconditioner: returns ``apply(r) = r / diag(A)`` as a callable."""
    d = [A[i][i] for i in range(len(A))]
    if any(abs(di) < 1e-300 for di in d):
        raise ValueError("zero on the diagonal; Jacobi preconditioner undefined")
    return lambda r: [r[i] / d[i] for i in range(len(r))]


def incomplete_cholesky(A):
    """IC(0) incomplete-Cholesky factor ``L`` of an SPD matrix ``A`` (same sparsity, no fill-in).

    Returns the lower-triangular ``L`` with ``L L^T ~ A`` on the nonzero pattern of ``A``. Entries
    where ``A[i][j] == 0`` are kept zero (the "no fill" rule). Use with :func:`ic_apply` /
    :func:`preconditioned_cg`.
    """
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            if A[i][j] == 0.0 and i != j:
                continue                        # preserve sparsity: no fill
            s = A[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                if s <= 0:
                    s = 1e-12                    # guard against loss of definiteness
                L[i][j] = math.sqrt(s)
            else:
                L[i][j] = s / L[j][j]
    return L


def ic_apply(L, r):
    """Apply the IC preconditioner: solve ``L L^T z = r`` by forward/back substitution."""
    n = len(L)
    y = [0.0] * n
    for i in range(n):
        y[i] = (r[i] - sum(L[i][k] * y[k] for k in range(i))) / L[i][i]
    z = [0.0] * n
    for i in range(n - 1, -1, -1):
        z[i] = (y[i] - sum(L[k][i] * z[k] for k in range(i + 1, n))) / L[i][i]
    return z


def preconditioned_cg(A, b, apply_minv=None, x0=None, tol=1e-10, max_iter=None):
    """Preconditioned conjugate gradient for a symmetric positive-definite ``A``.

    ``apply_minv`` is a callable ``r -> M^{-1} r`` (default: identity = plain CG). Returns a dict
    with ``x``, ``residual_norm``, ``n_iter`` and ``converged``. With a good preconditioner it
    converges in far fewer iterations than unpreconditioned CG.
    """
    n = len(b)
    if apply_minv is None:
        apply_minv = lambda r: list(r)
    if max_iter is None:
        max_iter = 10 * n
    x = [0.0] * n if x0 is None else [float(v) for v in x0]

    def mv(v):
        return [sum(A[i][j] * v[j] for j in range(n)) for i in range(n)]

    r = [b[i] - mv(x)[i] for i in range(n)]
    z = apply_minv(r)
    p = list(z)
    rz = _dot(r, z)
    bnorm = math.sqrt(_dot(b, b)) or 1.0

    n_iter = 0
    for it in range(max_iter):
        n_iter = it + 1
        if math.sqrt(_dot(r, r)) / bnorm <= tol:
            break
        Ap = mv(p)
        alpha = rz / _dot(p, Ap)
        x = [x[i] + alpha * p[i] for i in range(n)]
        r = [r[i] - alpha * Ap[i] for i in range(n)]
        z = apply_minv(r)
        rz_new = _dot(r, z)
        beta = rz_new / rz if rz != 0 else 0.0
        p = [z[i] + beta * p[i] for i in range(n)]
        rz = rz_new
    rnorm = math.sqrt(_dot(r, r))
    return {"x": x, "residual_norm": rnorm, "n_iter": n_iter,
            "converged": rnorm / bnorm <= tol}
