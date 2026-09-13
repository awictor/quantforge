"""Shooting method for two-point boundary-value problems.

A boundary-value problem fixes the solution at *both* ends -- ``y(a) = alpha``,
``y(b) = beta`` for ``y'' = f(t, y, y')`` -- rather than giving initial slope. The
shooting method turns it back into an initial-value problem: guess the initial slope
``y'(a) = s``, integrate to ``b`` with an ODE solver, and adjust ``s`` (by a root find on
the terminal residual ``y(b; s) - beta``) until the far boundary is hit. Reuses the
adaptive RK45 integrator and Brent root finder. Pure standard library.
"""

from .ode import rk45
from .rootfind import brent


def shooting_bvp(f, a, b, alpha, beta, s_lo, s_hi, tol=1e-8):
    """Solve ``y'' = f(t, y, y')`` with ``y(a)=alpha``, ``y(b)=beta`` by shooting.

    ``f(t, y, yp)`` returns ``y''``. ``s_lo, s_hi`` bracket the unknown initial slope
    ``y'(a)`` (the terminal residual must change sign across them). Returns a dict with
    the found initial ``slope``, and the solution ``ts`` / ``ys`` (state = ``[y, y']``)
    from the accepted RK45 steps.
    """
    def terminal_y(s):
        # Integrate the first-order system [y, y'] from a to b with y'(a)=s.
        _, ys = rk45(lambda t, u: [u[1], f(t, u[0], u[1])], a, [alpha, s], b, tol=tol)
        return ys[-1][0]

    res_lo = terminal_y(s_lo) - beta
    res_hi = terminal_y(s_hi) - beta
    if res_lo * res_hi > 0:
        raise ValueError("initial-slope bracket does not straddle the boundary; "
                         "widen s_lo/s_hi")

    slope = brent(lambda s: terminal_y(s) - beta, s_lo, s_hi, tol=tol)
    ts, ys = rk45(lambda t, u: [u[1], f(t, u[0], u[1])], a, [alpha, slope], b, tol=tol)
    return {"slope": slope, "ts": ts, "ys": ys}
