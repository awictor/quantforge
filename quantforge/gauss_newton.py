"""Gauss-Newton / Levenberg-Marquardt least squares with exact autodiff Jacobians.

For a nonlinear least-squares problem ``min sum_i r_i(p)^2``, the Gauss-Newton method
approximates the Hessian by ``J^T J`` (dropping the second-derivative term, which is small near
a good fit) and solves ``J^T J delta = -J^T r`` for the parameter step. Levenberg-Marquardt
interpolates between Gauss-Newton and gradient descent via a damping term ``lambda diag(J^T J)``
that is grown when a step fails and shrunk when it succeeds.

Unlike :func:`quantforge.levenberg_marquardt`, whose residual Jacobian is taken by finite
differences, this evaluates the *exact* Jacobian with reverse-mode autodiff -- the residual
function is written once with :class:`quantforge.reverse_ad.Var` arithmetic. Pure standard
library.
"""

import math

from .reverse_jacobian import reverse_jacobian
from .reverse_ad import Var
from .lu import lu_solve


def _sse(residual, p):
    r = [ri.value for ri in residual([Var(pi) for pi in p])]
    return sum(ri * ri for ri in r)


def gauss_newton(residual, p0, tol=1e-10, max_iter=200):
    """Minimize ``sum residual(p)^2`` by Levenberg-Marquardt with an exact autodiff Jacobian.

    ``residual`` takes a length-``k`` list of :class:`quantforge.reverse_ad.Var` and returns a
    list of ``m`` residual :class:`Var` (``m >= k``). Returns a dict with ``p`` (the fitted
    parameters), ``cost`` (half the sum of squared residuals ``0.5 * ||r||^2``), ``n_iter``,
    ``converged`` and ``grad_norm`` (norm of ``J^T r``).
    """
    k = len(p0)
    if k == 0:
        raise ValueError("need at least one parameter")
    p = [float(v) for v in p0]
    lam = 1e-3
    cost = 0.5 * _sse(residual, p)

    n_iter = 0
    for it in range(max_iter):
        n_iter = it + 1
        # exact Jacobian J (m x k) and residual vector r (m)
        inputs = [Var(pi) for pi in p]
        rvars = residual(inputs)
        r = [rv.value for rv in rvars]
        J = reverse_jacobian(residual, p)
        m = len(r)

        # gradient g = J^T r  and normal-equation matrix A = J^T J
        g = [sum(J[i][a] * r[i] for i in range(m)) for a in range(k)]
        gnorm = math.sqrt(sum(gi * gi for gi in g))
        if gnorm < tol:
            return {"p": p, "cost": cost, "n_iter": n_iter, "converged": True,
                    "grad_norm": gnorm}
        A = [[sum(J[i][a] * J[i][b] for i in range(m)) for b in range(k)] for a in range(k)]

        # LM: grow lambda until a step reduces the cost
        improved = False
        for _ in range(40):
            damped = [[A[a][b] + (lam * A[a][a] if a == b else 0.0) for b in range(k)]
                      for a in range(k)]
            try:
                delta = lu_solve(damped, [-gi for gi in g])
            except (ZeroDivisionError, ValueError):
                lam *= 10.0
                continue
            p_new = [p[a] + delta[a] for a in range(k)]
            cost_new = 0.5 * _sse(residual, p_new)
            if cost_new < cost:
                p, cost = p_new, cost_new
                lam = max(lam * 0.5, 1e-12)   # accept -> trust GN more next time
                improved = True
                break
            lam *= 10.0                        # reject -> more gradient-descent-like
        if not improved:
            # cannot reduce further; treat as converged at a local min
            return {"p": p, "cost": cost, "n_iter": n_iter, "converged": True,
                    "grad_norm": gnorm}

    inputs = [Var(pi) for pi in p]
    J = reverse_jacobian(residual, p)
    r = [rv.value for rv in residual(inputs)]
    m = len(r)
    g = [sum(J[i][a] * r[i] for i in range(m)) for a in range(k)]
    gnorm = math.sqrt(sum(gi * gi for gi in g))
    return {"p": p, "cost": cost, "n_iter": n_iter, "converged": gnorm < tol,
            "grad_norm": gnorm}
