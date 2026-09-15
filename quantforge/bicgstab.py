"""BiCGSTAB: short-recurrence Krylov solver for nonsymmetric systems.

GMRES (:func:`quantforge.gmres.gmres`) minimizes the residual exactly but stores the whole
Krylov basis, so its memory and per-step cost grow with the iteration count (hence restarting).
BiCGSTAB (van der Vorst 1992) is the popular *short-recurrence* alternative: it solves the same
general ``A x = b`` with a fixed, small amount of work and storage per step -- a few vectors, no
growing basis -- by combining the bi-conjugate-gradient recurrence with a stabilizing GMRES(1)
step. It often converges as fast as GMRES on well-conditioned nonsymmetric systems (discretized
convection-diffusion, Newton steps) at constant memory.

Uses ``A`` only through matrix-vector products (matrix or callable). Pure standard library.
"""

import math


def _as_matvec(A):
    if callable(A):
        return A
    return lambda v: [sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A))]


def _dot(a, b):
    return sum(a[i] * b[i] for i in range(len(a)))


def bicgstab(A, b, x0=None, tol=1e-10, max_iter=None):
    """Solve ``A x = b`` by BiCGSTAB for a general non-singular operator ``A``.

    ``A`` is a matrix or a callable ``v -> A@v``. Returns a dict with ``x``, ``residuals`` (the
    relative residual norm each iteration), ``converged`` and ``n_iter``. Constant memory: it
    keeps only a handful of length-``n`` vectors regardless of iteration count.
    """
    n = len(b)
    matvec = _as_matvec(A)
    x = [0.0] * n if x0 is None else [float(v) for v in x0]
    bnorm = math.sqrt(_dot(b, b)) or 1.0
    if max_iter is None:
        max_iter = 10 * n

    Ax = matvec(x)
    r = [b[i] - Ax[i] for i in range(n)]
    r_hat = list(r)                       # shadow residual (fixed)
    residuals = [math.sqrt(_dot(r, r)) / bnorm]
    if residuals[0] <= tol:
        return {"x": x, "residuals": residuals, "converged": True, "n_iter": 0}

    rho_old = alpha = omega = 1.0
    v = [0.0] * n
    p = [0.0] * n

    for it in range(1, max_iter + 1):
        rho = _dot(r_hat, r)
        if abs(rho) < 1e-300:
            break                          # breakdown
        beta = (rho / rho_old) * (alpha / omega)
        p = [r[i] + beta * (p[i] - omega * v[i]) for i in range(n)]
        v = matvec(p)
        denom = _dot(r_hat, v)
        if abs(denom) < 1e-300:
            break
        alpha = rho / denom
        s = [r[i] - alpha * v[i] for i in range(n)]
        s_norm = math.sqrt(_dot(s, s))
        if s_norm / bnorm <= tol:
            x = [x[i] + alpha * p[i] for i in range(n)]
            residuals.append(s_norm / bnorm)
            return {"x": x, "residuals": residuals, "converged": True, "n_iter": it}
        t = matvec(s)
        tt = _dot(t, t)
        omega = _dot(t, s) / tt if tt > 1e-300 else 0.0
        x = [x[i] + alpha * p[i] + omega * s[i] for i in range(n)]
        r = [s[i] - omega * t[i] for i in range(n)]
        res = math.sqrt(_dot(r, r)) / bnorm
        residuals.append(res)
        if res <= tol:
            return {"x": x, "residuals": residuals, "converged": True, "n_iter": it}
        if abs(omega) < 1e-300:
            break
        rho_old = rho

    return {"x": x, "residuals": residuals, "converged": residuals[-1] <= tol, "n_iter": len(residuals) - 1}
