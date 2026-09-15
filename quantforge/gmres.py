"""GMRES: Krylov solver for general (nonsymmetric) linear systems.

Conjugate gradient (:func:`quantforge.conjugate_gradient.conjugate_gradient`) solves ``A x = b``
only when ``A`` is symmetric positive-definite. GMRES (Saad & Schultz 1986) handles *any*
non-singular ``A``: it builds an orthonormal Krylov basis by Arnoldi iteration and, at each step,
picks the solution in that subspace minimizing the residual ``||b - A x||`` via a small
least-squares problem (solved incrementally with Givens rotations on the Hessenberg matrix). Like
CG it uses ``A`` only through matrix-vector products, so it suits large sparse / matrix-free
operators; it is the standard solver for the nonsymmetric systems from discretized PDEs and
Newton steps.

This is (restarted-capable) GMRES with the Arnoldi process and incremental Givens QR of the
Hessenberg matrix, returning the solution and residual history. Pure standard library.
"""

import math


def _as_matvec(A):
    if callable(A):
        return A
    return lambda v: [sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A))]


def _dot(a, b):
    return sum(a[i] * b[i] for i in range(len(a)))


def gmres(A, b, x0=None, tol=1e-10, max_iter=None, restart=None):
    """Solve ``A x = b`` by GMRES for a general non-singular operator ``A``.

    ``A`` is an ``n x n`` matrix or a callable ``v -> A@v``. ``restart`` sets the Krylov subspace
    size before restarting (default: ``min(n, max_iter or n)``). Returns a dict with ``x`` (the
    solution), ``residuals`` (the residual norm after each outer iteration), ``converged`` and
    ``n_iter`` (total inner iterations).
    """
    n = len(b)
    matvec = _as_matvec(A)
    x = [0.0] * n if x0 is None else [float(v) for v in x0]
    bnorm = math.sqrt(_dot(b, b)) or 1.0
    max_outer = max_iter if max_iter is not None else 10 * n
    m = restart if restart is not None else min(n, max_outer)

    residuals = []
    total_inner = 0
    for _outer in range(max_outer):
        # r0 = b - A x
        Ax = matvec(x)
        r = [b[i] - Ax[i] for i in range(n)]
        beta = math.sqrt(_dot(r, r))
        residuals.append(beta / bnorm)
        if beta / bnorm <= tol:
            return {"x": x, "residuals": residuals, "converged": True, "n_iter": total_inner}

        V = [[ri / beta for ri in r]]        # Arnoldi basis
        H = []                                # Hessenberg columns (each length up to m+1)
        cs = []                               # Givens cosines / sines
        sn = []
        g = [beta] + [0.0] * m                # RHS of the least-squares problem

        k = 0
        for k in range(m):
            total_inner += 1
            w = matvec(V[k])
            hcol = [0.0] * (k + 2)
            for i in range(k + 1):
                hcol[i] = _dot(w, V[i])
                w = [w[j] - hcol[i] * V[i][j] for j in range(n)]
            hcol[k + 1] = math.sqrt(_dot(w, w))
            if hcol[k + 1] > 1e-14:
                V.append([wj / hcol[k + 1] for wj in w])
            # apply previous Givens rotations to the new column
            for i in range(k):
                t = cs[i] * hcol[i] + sn[i] * hcol[i + 1]
                hcol[i + 1] = -sn[i] * hcol[i] + cs[i] * hcol[i + 1]
                hcol[i] = t
            # new Givens rotation to zero hcol[k+1]
            denom = math.hypot(hcol[k], hcol[k + 1])
            c = hcol[k] / denom
            s = hcol[k + 1] / denom
            cs.append(c)
            sn.append(s)
            hcol[k] = c * hcol[k] + s * hcol[k + 1]
            hcol[k + 1] = 0.0
            H.append(hcol)
            # rotate the RHS
            g[k + 1] = -s * g[k]
            g[k] = c * g[k]
            residuals.append(abs(g[k + 1]) / bnorm)
            if abs(g[k + 1]) / bnorm <= tol or hcol[k + 1] == 0.0 and False:
                k += 1
                break
            if abs(g[k + 1]) / bnorm <= tol:
                k += 1
                break
        # solve the k x k upper-triangular system H y = g and update x
        kk = k
        y = [0.0] * kk
        for i in range(kk - 1, -1, -1):
            s = g[i] - sum(H[j][i] * y[j] for j in range(i + 1, kk))
            y[i] = s / H[i][i]
        for i in range(n):
            x[i] += sum(y[j] * V[j][i] for j in range(kk))

        Ax = matvec(x)
        r = [b[i] - Ax[i] for i in range(n)]
        if math.sqrt(_dot(r, r)) / bnorm <= tol:
            residuals.append(math.sqrt(_dot(r, r)) / bnorm)
            return {"x": x, "residuals": residuals, "converged": True, "n_iter": total_inner}

    Ax = matvec(x)
    r = [b[i] - Ax[i] for i in range(n)]
    resid = math.sqrt(_dot(r, r)) / bnorm
    return {"x": x, "residuals": residuals, "converged": resid <= tol, "n_iter": total_inner}
