"""Multivariate autodiff helpers: exact gradients and derivative-free Newton.

Built on the forward-mode :class:`quantforge.dual.Dual` numbers:

  * ``dual_gradient`` -- the exact gradient of a scalar function of a vector, computed
    by seeding one component's derivative to 1 at a time (``n`` evaluations, each
    exact -- no finite-difference error).
  * ``dual_newton`` -- Newton's method for a scalar root that needs *only* the function:
    the derivative at each step comes from autodiff, so there is no hand-coded or
    finite-difference ``f'`` to supply or get wrong.

Pure standard library.
"""

from .dual import Dual


def dual_gradient(f, x):
    """Exact gradient of a scalar ``f`` of a vector ``x`` by forward-mode autodiff.

    ``f`` takes a list of (possibly :class:`Dual`) components and returns a scalar
    ``Dual``. Seeds each coordinate's derivative to 1 in turn, so the result is the
    vector of partial derivatives with no truncation error. ``n`` function evaluations.
    """
    n = len(x)
    grad = []
    for j in range(n):
        args = [Dual(x[i], 1.0 if i == j else 0.0) for i in range(n)]
        out = f(args)
        grad.append(out.deriv if isinstance(out, Dual) else 0.0)
    return grad


def dual_newton(f, x0, tol=1e-12, max_iter=100):
    """Newton's method for ``f(x) = 0`` using autodiff for the derivative.

    ``f`` must accept a :class:`Dual` and return a :class:`Dual`. Each step evaluates
    ``f`` once on a seeded dual, reading both ``f(x)`` and ``f'(x)`` at no extra cost,
    so no separate derivative function is needed. Returns a dict with ``root``,
    ``iterations`` and ``converged``. Raises if the derivative vanishes.
    """
    x = float(x0)
    for it in range(1, max_iter + 1):
        d = f(Dual(x, 1.0))
        fx, fpx = d.value, d.deriv
        if fpx == 0.0:
            raise ValueError("zero derivative in Newton iteration")
        step = fx / fpx
        x -= step
        if abs(step) < tol:
            return {"root": x, "iterations": it, "converged": True}
    return {"root": x, "iterations": max_iter, "converged": False}
