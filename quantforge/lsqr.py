"""LSQR: Krylov least-squares solver for rectangular systems.

Conjugate gradient / GMRES / BiCGSTAB solve *square* systems. LSQR (Paige & Saunders 1982)
solves the least-squares problem ``min ||A x - b||`` for a general ``m x n`` matrix -- including
overdetermined (``m > n``, fitting) and underdetermined systems -- and is analytically
equivalent to conjugate gradient on the normal equations ``A^T A x = A^T b`` but far more
numerically stable, because it never forms ``A^T A`` (whose condition number is squared). It
works through the Golub-Kahan bidiagonalization, needing only products with ``A`` and ``A^T``.

This implements standard LSQR (unregularized), taking ``A`` as a matrix (its transpose is formed
internally) and returning the least-squares solution and residual norm. Pure standard library.
"""

import math


def _matvec(A, v):
    return [sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A))]


def _rmatvec(A, v):
    # A^T @ v
    n = len(A[0])
    return [sum(A[i][j] * v[i] for i in range(len(A))) for j in range(n)]


def _norm(v):
    return math.sqrt(sum(x * x for x in v))


def lsqr(A, b, tol=1e-10, max_iter=None):
    """Least-squares solution of ``min ||A x - b||`` by LSQR for an ``m x n`` matrix ``A``.

    Handles overdetermined (``m > n``) and square systems. Returns a dict with ``x`` (length
    ``n``), ``residual`` (``||A x - b||``), ``n_iter`` and ``converged`` (relative normal-equation
    residual ``||A^T r|| / (||A|| ||r||)`` below ``tol``).
    """
    m = len(A)
    n = len(A[0])
    if max_iter is None:
        max_iter = 4 * n

    # Golub-Kahan bidiagonalization initialization
    beta = _norm(b)
    if beta == 0.0:
        return {"x": [0.0] * n, "residual": 0.0, "n_iter": 0, "converged": True}
    u = [bi / beta for bi in b]
    v = _rmatvec(A, u)
    alpha = _norm(v)
    if alpha == 0.0:
        return {"x": [0.0] * n, "residual": beta, "n_iter": 0, "converged": True}
    v = [vi / alpha for vi in v]

    w = list(v)
    x = [0.0] * n
    phibar = beta
    rhobar = alpha
    anorm2 = alpha * alpha

    n_iter = 0
    for it in range(max_iter):
        n_iter = it + 1
        # continue bidiagonalization
        u2 = [_matvec(A, v)[i] - alpha * u[i] for i in range(m)]
        beta = _norm(u2)
        if beta > 0:
            u = [x_ / beta for x_ in u2]
            anorm2 += beta * beta
            v2 = [_rmatvec(A, u)[j] - beta * v[j] for j in range(n)]
            alpha = _norm(v2)
            if alpha > 0:
                v = [x_ / alpha for x_ in v2]
                anorm2 += alpha * alpha
        # orthogonal transformation (Givens) on the bidiagonal
        rho = math.hypot(rhobar, beta)
        c = rhobar / rho
        s = beta / rho
        theta = s * alpha
        rhobar = -c * alpha
        phi = c * phibar
        phibar = s * phibar
        # update x and the search direction w
        t1 = phi / rho
        t2 = -theta / rho
        x = [x[j] + t1 * w[j] for j in range(n)]
        w = [v[j] + t2 * w[j] for j in range(n)]
        # stopping test: ||A^T r|| = |phibar * alpha * c| ; ||r|| = phibar
        atr = abs(phibar * alpha * c)
        anorm = math.sqrt(anorm2)
        rnorm = phibar
        if rnorm < 1e-300 or atr / (anorm * rnorm + 1e-300) <= tol:
            break

    # actual residual
    r = [b[i] - _matvec(A, x)[i] for i in range(m)]
    rnorm = _norm(r)
    atr = _norm(_rmatvec(A, r))
    anorm = math.sqrt(anorm2)
    converged = atr / (anorm * rnorm + 1e-300) <= tol if rnorm > 0 else True
    return {"x": x, "residual": rnorm, "n_iter": n_iter, "converged": converged}
