"""Jacobians and Hessians built on reverse-mode autodiff.

:mod:`quantforge.reverse_ad` gives an exact gradient of a scalar function in one backward pass.
This module composes that primitive into the two matrices most numerical code needs:

* :func:`reverse_jacobian` -- for a vector-valued ``f: R^n -> R^m``, one backward pass per
  output component yields the exact ``m x n`` Jacobian (no differencing error).
* :func:`reverse_hessian` -- for a scalar ``f: R^n -> R``, the exact reverse gradient is
  differenced once by central finite differences to give the ``n x n`` Hessian. Because only
  the *second* derivative carries truncation error (the gradient is exact), this is far more
  accurate than differencing the function value twice, and the result is symmetrized.

These complement the fully numerical `jacobian`/`hessian` in :mod:`quantforge.numdiff`. Pure
standard library.
"""

from .reverse_ad import Var


def reverse_jacobian(f, x):
    """Exact Jacobian of ``f: R^n -> R^m`` at ``x`` via reverse-mode autodiff.

    ``f`` takes a list of :class:`Var` and returns a list of :class:`Var` (length ``m``).
    Returns the ``m x n`` Jacobian as a list of rows; row ``i`` is ``grad f_i``. Each output
    component is differentiated by its own backward pass.
    """
    n = len(x)
    inputs = [Var(xi) for xi in x]
    outputs = f(inputs)
    if isinstance(outputs, Var):
        outputs = [outputs]
    jac = []
    for out in outputs:
        # backward() only resets grads on the subgraph reachable from `out`; an input that
        # does not feed this output would otherwise keep a stale grad from a prior pass.
        for inp in inputs:
            inp.grad = 0.0
        out.backward()
        jac.append([inputs[j].grad for j in range(n)])
    return jac


def reverse_gradient_vector(f, x):
    """Exact gradient of scalar ``f`` at ``x`` (thin wrapper returning a plain list).

    ``f`` takes a list of :class:`Var` and returns a single :class:`Var`.
    """
    inputs = [Var(xi) for xi in x]
    out = f(inputs)
    out.backward()
    return [v.grad for v in inputs]


def reverse_hessian(f, x, h=1e-5):
    """Hessian of scalar ``f: R^n -> R`` at ``x`` by differencing the exact reverse gradient.

    ``f`` takes a list of :class:`Var` and returns a single :class:`Var`. Each column is a
    central difference of the (exact) gradient, ``(grad f(x + h e_j) - grad f(x - h e_j)) /
    (2 h)``; the result is symmetrized. Only the outer difference carries truncation error, so
    accuracy is far better than second-differencing ``f`` itself.
    """
    n = len(x)
    hess = [[0.0] * n for _ in range(n)]
    for j in range(n):
        xp = list(x)
        xm = list(x)
        xp[j] += h
        xm[j] -= h
        gp = reverse_gradient_vector(f, xp)
        gm = reverse_gradient_vector(f, xm)
        for i in range(n):
            hess[i][j] = (gp[i] - gm[i]) / (2.0 * h)
    # symmetrize (averages out the O(h^2) asymmetry)
    for i in range(n):
        for j in range(i + 1, n):
            avg = 0.5 * (hess[i][j] + hess[j][i])
            hess[i][j] = hess[j][i] = avg
    return hess
