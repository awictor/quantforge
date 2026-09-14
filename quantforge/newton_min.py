"""Newton's method for unconstrained minimization (with Hessian damping).

Gradient methods like :func:`quantforge.lbfgs.lbfgs` use only first derivatives. When the
Hessian is available -- and here it comes exactly from reverse-mode autodiff -- Newton's method
converges *quadratically* near a minimum: it steps to the minimum of the local quadratic model
``x <- x - H^{-1} g``. Far from the optimum the raw Newton step can point uphill (indefinite
Hessian), so this adds Levenberg-style diagonal damping ``(H + lambda I)`` until the step is a
descent direction, plus an Armijo backtracking line search for global robustness.

The objective is written with :class:`quantforge.reverse_ad.Var` arithmetic, so the gradient
and Hessian are obtained without any hand-coded derivatives. Pure standard library.
"""

import math

from .reverse_jacobian import reverse_gradient_vector, reverse_hessian
from .lu import lu_solve


def _dot(a, b):
    return sum(ai * bi for ai, bi in zip(a, b))


def newton_min(func, x0, tol=1e-9, max_iter=100, h=1e-5):
    """Minimize ``func`` from ``x0`` by damped Newton with exact autodiff derivatives.

    ``func`` takes a list of :class:`quantforge.reverse_ad.Var` and returns a single ``Var``.
    The gradient is exact (reverse mode); the Hessian is the reverse gradient differenced once
    (:func:`quantforge.reverse_jacobian.reverse_hessian`). Returns a dict with ``x``, ``fun``,
    ``n_iter``, ``converged`` (gradient norm below ``tol``) and ``grad_norm``.
    """
    n = len(x0)
    if n == 0:
        raise ValueError("need at least one dimension")
    x = [float(v) for v in x0]

    def fval(z):
        from .reverse_ad import Var
        return func([Var(zi) for zi in z]).value

    fx = fval(x)
    n_iter = 0
    for it in range(max_iter):
        n_iter = it + 1
        g = reverse_gradient_vector(func, x)
        gnorm = math.sqrt(_dot(g, g))
        if gnorm < tol:
            return {"x": x, "fun": fx, "n_iter": n_iter, "converged": True,
                    "grad_norm": gnorm}
        H = reverse_hessian(func, x, h)

        # damped Newton: grow lambda until (H + lambda I) gives a descent direction
        lam = 0.0
        step = None
        for _ in range(40):
            damped = [[H[i][j] + (lam if i == j else 0.0) for j in range(n)] for i in range(n)]
            try:
                p = lu_solve(damped, [-gi for gi in g])
            except (ZeroDivisionError, ValueError):
                lam = max(1e-3, lam * 10.0)
                continue
            if _dot(g, p) < 0:      # descent direction found
                step = p
                break
            lam = max(1e-3, lam * 10.0)
        if step is None:
            step = [-gi for gi in g]   # fall back to steepest descent

        # Armijo backtracking line search along the (damped) Newton direction
        slope = _dot(g, step)
        alpha = 1.0
        c1 = 1e-4
        while True:
            x_new = [x[i] + alpha * step[i] for i in range(n)]
            f_new = fval(x_new)
            if f_new <= fx + c1 * alpha * slope or alpha < 1e-14:
                break
            alpha *= 0.5
        x, fx = x_new, f_new

    g = reverse_gradient_vector(func, x)
    gnorm = math.sqrt(_dot(g, g))
    return {"x": x, "fun": fx, "n_iter": n_iter, "converged": gnorm < tol,
            "grad_norm": gnorm}
