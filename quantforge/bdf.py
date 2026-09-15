"""BDF2: second-order backward-differentiation formula for stiff ODEs.

Backward Euler (:func:`quantforge.stiff_ode.backward_euler`) is A-stable but only first order, so
it is heavily damped and inaccurate. BDF2 is the second-order backward-differentiation formula --
the workhorse multistep method for stiff systems (the core of production stiff solvers). It fits
a quadratic through the last two states and the new one:

    (3 y_{n+1} - 4 y_n + y_{n-1}) / (2 h) = f(t_{n+1}, y_{n+1}) .

BDF2 is A(alpha)-stable with alpha ~ 90 degrees (effectively A-stable for practical problems) and
second order, so it steps by accuracy on stiff problems while damping the fast transient. The
first step is bootstrapped with backward Euler; each step solves the implicit equation by
Newton's method with a finite-difference Jacobian and :func:`quantforge.lu.lu_solve`. Pure
standard library.
"""

from .lu import lu_solve
from .stiff_ode import _as_vec, _fd_jacobian, _newton_step, backward_euler


def bdf2(f, y0, t0, t1, steps):
    """Integrate ``y' = f(t, y)`` from ``t0`` to ``t1`` in ``steps`` BDF2 steps.

    ``f(t, y)`` returns the derivative (scalar or list); ``y0`` matches. Second order, A(alpha)-
    stable. The first step uses backward Euler to seed the two-step history. Returns ``(ts, ys)``
    with states as lists.
    """
    if steps < 2:
        raise ValueError("BDF2 needs at least 2 steps")
    scalar = not isinstance(y0, (list, tuple))
    h = (t1 - t0) / steps
    n = len(_as_vec(y0))

    def f_vec(t, yv):
        r = f(t, yv[0]) if scalar else f(t, yv)
        return _as_vec(r)

    # bootstrap: one backward-Euler step to get y_1
    _, be = backward_euler(f, y0, t0, t0 + h, 1)
    ys = [list(be[0]), list(be[1])]
    ts = [t0, t0 + h]

    for k in range(1, steps):
        t_next = t0 + (k + 1) * h
        y_nm1 = ys[-2]
        y_n = ys[-1]

        def residual(yn):
            fn = f_vec(t_next, yn)
            return [3.0 * yn[i] - 4.0 * y_n[i] + y_nm1[i] - 2.0 * h * fn[i] for i in range(n)]

        def jac_of_res(yn):
            J = _fd_jacobian(f_vec, t_next, yn)
            return [[(3.0 if i == j else 0.0) - 2.0 * h * J[i][j] for j in range(n)]
                    for i in range(n)]

        y = _newton_step(residual, jac_of_res, y_n)
        ys.append(list(y))
        ts.append(t_next)
    return ts, ys
