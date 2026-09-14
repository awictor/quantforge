"""Find all roots of a function on an interval by sign-change scanning.

A single bracketing solver (Brent) needs a sign change to start; to find *every* root of
a nonlinear function on ``[a, b]`` you first scan a fine grid for sign changes, then
refine each bracket. This locates all simple (sign-changing) roots -- the practical tool
for solving a transcendental equation, finding implied levels, or seeding a global
optimizer. Even-multiplicity roots (no sign change) are missed; use a denser grid or add
the function's derivative if those matter. Pure standard library.
"""

from .rootfind import brent


def find_all_roots(f, a, b, n=1000, tol=1e-12):
    """All sign-changing roots of ``f`` on ``[a, b]`` via a grid scan plus Brent refinement.

    Splits ``[a, b]`` into ``n`` subintervals, and wherever ``f`` changes sign (or hits
    exactly zero at a node) brackets and refines a root with Brent's method. Returns the
    roots in increasing order, de-duplicated. Increase ``n`` to catch roots closer than
    the grid spacing.
    """
    if b <= a:
        raise ValueError("require a < b")
    if n < 1:
        raise ValueError("n must be >= 1")
    step = (b - a) / n
    roots = []
    x_prev = a
    f_prev = f(x_prev)
    if f_prev == 0.0:
        roots.append(x_prev)
    for i in range(1, n + 1):
        x_cur = a + i * step
        f_cur = f(x_cur)
        if f_cur == 0.0:
            roots.append(x_cur)
        elif f_prev * f_cur < 0.0:
            roots.append(brent(f, x_prev, x_cur, tol=tol))
        x_prev, f_prev = x_cur, f_cur
    # De-duplicate roots that landed within tolerance of one another.
    roots.sort()
    unique = []
    for r in roots:
        if not unique or abs(r - unique[-1]) > max(tol * 10, 1e-9):
            unique.append(r)
    return unique


def count_sign_changes(f, a, b, n=1000):
    """Number of sign changes of ``f`` on ``[a, b]`` over an ``n``-point grid.

    A quick lower bound on the number of simple roots (each sign change brackets at least
    one). Cheaper than :func:`find_all_roots` when only the count is needed.
    """
    if b <= a:
        raise ValueError("require a < b")
    if n < 1:
        raise ValueError("n must be >= 1")
    step = (b - a) / n
    changes = 0
    f_prev = f(a)
    for i in range(1, n + 1):
        f_cur = f(a + i * step)
        if f_prev * f_cur < 0.0:
            changes += 1
        f_prev = f_cur
    return changes
