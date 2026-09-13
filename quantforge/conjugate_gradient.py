"""Iterative linear solvers: conjugate gradient, Gauss-Seidel, Jacobi.

Direct methods (LU, Cholesky) cost ``O(n^3)``; for large sparse or structured systems
an iterative solver that only needs matrix-vector products is far cheaper:

  * ``conjugate_gradient`` -- the method of choice for a symmetric positive-definite
    ``A``. In exact arithmetic it converges in at most ``n`` steps and usually far
    fewer, driving the residual down along ``A``-conjugate directions.
  * ``gauss_seidel`` / ``jacobi`` -- classic stationary iterations that converge for
    diagonally dominant systems; Gauss-Seidel uses each freshly updated component
    immediately, so it typically converges about twice as fast as Jacobi.

Each returns the solution, the residual norm, and the iteration count. Pure standard
library.
"""

import math


def _matvec(A, x):
    return [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]


def conjugate_gradient(A, b, x0=None, tol=1e-10, max_iter=None):
    """Solve a symmetric positive-definite system ``A x = b`` by conjugate gradient.

    Returns a dict with ``x`` (solution), ``residual_norm`` (``||b - A x||``) and
    ``n_iter``. ``A`` must be symmetric positive-definite for convergence; converges in
    at most ``n`` iterations in exact arithmetic.
    """
    n = len(b)
    if len(A) != n:
        raise ValueError("A and b dimensions must match")
    if max_iter is None:
        max_iter = 10 * n
    x = [0.0] * n if x0 is None else list(x0)
    r = [b[i] - sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]
    p = list(r)
    rs_old = sum(ri * ri for ri in r)
    if math.sqrt(rs_old) < tol:
        return {"x": x, "residual_norm": math.sqrt(rs_old), "n_iter": 0}

    n_iter = 0
    for it in range(1, max_iter + 1):
        n_iter = it
        Ap = _matvec(A, p)
        pAp = sum(p[i] * Ap[i] for i in range(n))
        if pAp == 0.0:
            break
        alpha = rs_old / pAp
        x = [x[i] + alpha * p[i] for i in range(n)]
        r = [r[i] - alpha * Ap[i] for i in range(n)]
        rs_new = sum(ri * ri for ri in r)
        if math.sqrt(rs_new) < tol:
            rs_old = rs_new
            break
        beta = rs_new / rs_old
        p = [r[i] + beta * p[i] for i in range(n)]
        rs_old = rs_new
    return {"x": x, "residual_norm": math.sqrt(rs_old), "n_iter": n_iter}


def _stationary(A, b, use_latest, tol, max_iter):
    n = len(b)
    if len(A) != n:
        raise ValueError("A and b dimensions must match")
    if any(A[i][i] == 0.0 for i in range(n)):
        raise ValueError("zero on the diagonal; reorder the system")
    if max_iter is None:
        max_iter = 1000
    x = [0.0] * n
    n_iter = 0
    for it in range(1, max_iter + 1):
        n_iter = it
        x_new = list(x)
        for i in range(n):
            src = x_new if use_latest else x
            s = sum(A[i][j] * src[j] for j in range(n) if j != i)
            x_new[i] = (b[i] - s) / A[i][i]
        # Residual on the updated iterate.
        resid = [b[i] - sum(A[i][j] * x_new[j] for j in range(n)) for i in range(n)]
        rn = math.sqrt(sum(rr * rr for rr in resid))
        x = x_new
        if rn < tol:
            return {"x": x, "residual_norm": rn, "n_iter": n_iter}
    resid = [b[i] - sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]
    return {"x": x, "residual_norm": math.sqrt(sum(rr * rr for rr in resid)),
            "n_iter": n_iter}


def gauss_seidel(A, b, tol=1e-10, max_iter=None):
    """Gauss-Seidel iteration for ``A x = b`` (uses freshly-updated components).

    Converges for diagonally dominant or SPD ``A``. Returns ``x``, ``residual_norm``
    and ``n_iter``.
    """
    return _stationary(A, b, use_latest=True, tol=tol, max_iter=max_iter)


def jacobi(A, b, tol=1e-10, max_iter=None):
    """Jacobi iteration for ``A x = b`` (all components updated from the old iterate).

    Converges for diagonally dominant ``A``; slower than Gauss-Seidel. Returns ``x``,
    ``residual_norm`` and ``n_iter``.
    """
    return _stationary(A, b, use_latest=False, tol=tol, max_iter=max_iter)
