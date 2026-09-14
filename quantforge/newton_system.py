"""Root finding for systems of nonlinear equations: Newton and Broyden.

Where :mod:`quantforge.rootfind` solves one scalar equation, these solve a vector system
``F(x) = 0`` in ``n`` unknowns. ``newton_system`` uses the finite-difference Jacobian and
a linear solve each step (quadratic convergence, but ``n`` extra function evaluations per
iteration for the Jacobian). ``broyden`` maintains a quasi-Newton approximation of the
inverse Jacobian, updating it rank-1 from each step -- cheaper per iteration when
function evaluations dominate. Pure standard library.
"""

from .numdiff import jacobian
from .lu import lu_solve


def _norm(v):
    return sum(c * c for c in v) ** 0.5


def newton_system(f, x0, tol=1e-10, max_iter=100, rel_step=1e-6):
    """Solve ``f(x) = 0`` for a vector function by Newton's method.

    ``f`` maps a length-``n`` list to a length-``n`` list; ``x0`` is the initial guess.
    Each step solves ``J dx = -f(x)`` with the central-difference Jacobian ``J`` and
    updates ``x += dx``. Returns ``(solution, iterations)``. Raises if the Jacobian is
    singular or convergence is not reached within ``max_iter``.
    """
    x = [float(v) for v in x0]
    for it in range(1, max_iter + 1):
        fx = f(x)
        if _norm(fx) < tol:
            return x, it - 1
        J = jacobian(f, x, rel_step=rel_step)
        neg = [-v for v in fx]
        try:
            dx = lu_solve(J, neg)
        except (ValueError, ZeroDivisionError):
            raise ValueError("singular Jacobian at iteration %d" % it)
        x = [x[i] + dx[i] for i in range(len(x))]
        if _norm(dx) < tol:
            return x, it
    raise ValueError("Newton did not converge in %d iterations" % max_iter)


def broyden(f, x0, tol=1e-10, max_iter=200, rel_step=1e-6):
    """Solve ``f(x) = 0`` by Broyden's (good) quasi-Newton method.

    Seeds the inverse-Jacobian estimate from one finite-difference Jacobian, then updates
    it rank-1 from each step's secant equation -- avoiding a fresh Jacobian per iteration.
    Returns ``(solution, iterations)``. Raises if the seed Jacobian is singular or it
    fails to converge.
    """
    x = [float(v) for v in x0]
    n = len(x)
    fx = f(x)
    if _norm(fx) < tol:
        return x, 0
    # Seed inverse Jacobian H = J^{-1} via solving J H = I column by column.
    J = jacobian(f, x, rel_step=rel_step)
    H = _inverse(J)
    for it in range(1, max_iter + 1):
        # step dx = -H fx
        dx = [-sum(H[i][k] * fx[k] for k in range(n)) for i in range(n)]
        x_new = [x[i] + dx[i] for i in range(n)]
        fx_new = f(x_new)
        if _norm(fx_new) < tol or _norm(dx) < tol:
            return x_new, it
        # Broyden good update of H.
        df = [fx_new[k] - fx[k] for k in range(n)]
        Hdf = [sum(H[i][k] * df[k] for k in range(n)) for i in range(n)]
        denom = sum(dx[i] * Hdf[i] for i in range(n))
        if denom == 0.0:
            raise ValueError("Broyden update broke down at iteration %d" % it)
        # H += ((dx - H df) / (dx^T H df)) dx^T H
        correction = [dx[i] - Hdf[i] for i in range(n)]
        dxH = [sum(dx[k] * H[k][j] for k in range(n)) for j in range(n)]
        for i in range(n):
            for j in range(n):
                H[i][j] += correction[i] * dxH[j] / denom
        x, fx = x_new, fx_new
    raise ValueError("Broyden did not converge in %d iterations" % max_iter)


def _inverse(A):
    """Inverse of a square float matrix via column solves (raises if singular)."""
    n = len(A)
    cols = []
    for c in range(n):
        e = [1.0 if i == c else 0.0 for i in range(n)]
        try:
            cols.append(lu_solve(A, e))
        except (ValueError, ZeroDivisionError):
            raise ValueError("singular Jacobian")
    return [[cols[c][r] for c in range(n)] for r in range(n)]
