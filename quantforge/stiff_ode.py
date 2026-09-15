"""Implicit ODE integrators for stiff systems (backward Euler, trapezoidal).

Explicit methods (:mod:`quantforge.ode` RK4/RK45, :mod:`quantforge.multistep_ode`) must take
tiny steps on *stiff* problems -- systems with widely separated time scales, like a fast decay
coupled to slow dynamics -- or they blow up. Implicit methods solve for the next state rather
than extrapolating to it, giving unconditional (A-)stability so the step size is set by accuracy,
not stability.

* :func:`backward_euler` -- ``y_{n+1} = y_n + h f(t_{n+1}, y_{n+1})``: A-stable, first order.
* :func:`trapezoidal` -- ``y_{n+1} = y_n + h/2 (f_n + f_{n+1})``: A-stable, second order.

Each step solves the implicit equation by Newton's method with a finite-difference Jacobian and
:func:`quantforge.lu.lu_solve`. State is a scalar or a vector (list). Pure standard library.
"""

from .lu import lu_solve


def _as_vec(y):
    return [y] if not isinstance(y, (list, tuple)) else list(y)


def _fd_jacobian(f, t, y, eps=1e-7):
    n = len(y)
    f0 = f(t, y)
    J = [[0.0] * n for _ in range(n)]
    for j in range(n):
        yp = list(y)
        h = eps * (abs(y[j]) + eps)
        yp[j] += h
        fj = f(t, yp)
        for i in range(n):
            J[i][j] = (fj[i] - f0[i]) / h
    return J


def _newton_step(residual, jac_of_res, y_guess, tol=1e-12, max_iter=50):
    y = list(y_guess)
    n = len(y)
    for _ in range(max_iter):
        r = residual(y)
        if max(abs(v) for v in r) < tol:
            break
        J = jac_of_res(y)
        delta = lu_solve(J, [-ri for ri in r])
        y = [y[i] + delta[i] for i in range(n)]
    return y


def backward_euler(f, y0, t0, t1, steps):
    """Integrate ``y' = f(t, y)`` from ``t0`` to ``t1`` in ``steps`` backward-Euler steps.

    ``f(t, y)`` returns the derivative (scalar or list); ``y0`` matches. A-stable, first order.
    Returns ``(ts, ys)`` -- the time points and states (states are lists even for scalar input).
    """
    scalar = not isinstance(y0, (list, tuple))
    y = _as_vec(y0)
    n = len(y)
    h = (t1 - t0) / steps
    ts = [t0]
    ys = [list(y)]

    def f_vec(t, yv):
        r = f(t, yv[0]) if scalar else f(t, yv)
        return _as_vec(r)

    for k in range(steps):
        t_next = t0 + (k + 1) * h
        y_cur = list(y)

        def residual(yn):
            fn = f_vec(t_next, yn)
            return [yn[i] - y_cur[i] - h * fn[i] for i in range(n)]

        def jac_of_res(yn):
            J = _fd_jacobian(f_vec, t_next, yn)
            return [[(1.0 if i == j else 0.0) - h * J[i][j] for j in range(n)] for i in range(n)]

        y = _newton_step(residual, jac_of_res, y_cur)
        ts.append(t_next)
        ys.append(list(y))
    return ts, ys


def trapezoidal(f, y0, t0, t1, steps):
    """Integrate ``y' = f(t, y)`` by the implicit trapezoidal rule (A-stable, second order).

    ``y_{n+1} = y_n + h/2 (f(t_n, y_n) + f(t_{n+1}, y_{n+1}))``. Same interface as
    :func:`backward_euler`.
    """
    scalar = not isinstance(y0, (list, tuple))
    y = _as_vec(y0)
    n = len(y)
    h = (t1 - t0) / steps
    ts = [t0]
    ys = [list(y)]

    def f_vec(t, yv):
        r = f(t, yv[0]) if scalar else f(t, yv)
        return _as_vec(r)

    for k in range(steps):
        t_cur = t0 + k * h
        t_next = t0 + (k + 1) * h
        y_cur = list(y)
        f_cur = f_vec(t_cur, y_cur)

        def residual(yn):
            fn = f_vec(t_next, yn)
            return [yn[i] - y_cur[i] - 0.5 * h * (f_cur[i] + fn[i]) for i in range(n)]

        def jac_of_res(yn):
            J = _fd_jacobian(f_vec, t_next, yn)
            return [[(1.0 if i == j else 0.0) - 0.5 * h * J[i][j] for j in range(n)]
                    for i in range(n)]

        y = _newton_step(residual, jac_of_res, y_cur)
        ts.append(t_next)
        ys.append(list(y))
    return ts, ys
