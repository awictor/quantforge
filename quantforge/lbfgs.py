"""Limited-memory BFGS (L-BFGS) unconstrained minimizer.

BFGS (:mod:`quantforge.bfgs`) stores the full ``n x n`` inverse-Hessian approximation, which
costs ``O(n^2)`` memory -- fine for tens of variables, painful for thousands. L-BFGS keeps only
the last ``m`` correction pairs ``(s, y)`` and reconstructs the search direction with Nocedal's
two-loop recursion in ``O(m n)`` time and memory, making it the standard choice for
large-scale smooth optimization (logistic regression, CRFs, calibration).

The line search enforces the strong Wolfe conditions so the curvature pairs stay
positive-definite. An analytic gradient may be supplied; otherwise central differences are
used. Pure standard library.
"""

import math


def _fd_grad(f, x, h=1e-6):
    n = len(x)
    g = [0.0] * n
    for i in range(n):
        xp = list(x)
        xm = list(x)
        xp[i] += h
        xm[i] -= h
        g[i] = (f(xp) - f(xm)) / (2.0 * h)
    return g


def _dot(a, b):
    return sum(ai * bi for ai, bi in zip(a, b))


def lbfgs(func, x0, grad=None, m=10, tol=1e-8, max_iter=500):
    """Minimize ``func`` from ``x0`` by L-BFGS with a strong-Wolfe line search.

    ``func`` maps a length-``n`` list to a scalar. ``grad`` is an optional gradient function
    (same signature, returning a length-``n`` list); if omitted, central differences are used.
    ``m`` is the history size. Returns a dict with ``x`` (minimizer), ``fun``, ``n_iter``,
    ``converged`` (gradient norm below ``tol``) and ``grad_norm``.
    """
    n = len(x0)
    if n == 0:
        raise ValueError("need at least one dimension")
    gfun = grad if grad is not None else (lambda z: _fd_grad(func, z))

    x = [float(v) for v in x0]
    fx = func(x)
    g = gfun(x)

    s_hist = []   # x_{k+1} - x_k
    y_hist = []   # g_{k+1} - g_k
    rho_hist = []

    n_iter = 0
    for it in range(max_iter):
        n_iter = it + 1
        gnorm = math.sqrt(_dot(g, g))
        if gnorm < tol:
            return {"x": x, "fun": fx, "n_iter": n_iter, "converged": True,
                    "grad_norm": gnorm}

        # ---- two-loop recursion: direction p = -H_k g ----
        q = list(g)
        alphas = []
        for s_i, y_i, rho_i in zip(reversed(s_hist), reversed(y_hist), reversed(rho_hist)):
            a = rho_i * _dot(s_i, q)
            alphas.append(a)
            q = [qj - a * yj for qj, yj in zip(q, y_i)]
        # initial Hessian scaling gamma = s.y / y.y
        if y_hist:
            sy = _dot(s_hist[-1], y_hist[-1])
            yy = _dot(y_hist[-1], y_hist[-1])
            gamma = sy / yy if yy > 0 else 1.0
        else:
            gamma = 1.0
        r = [gamma * qj for qj in q]
        for s_i, y_i, rho_i, a in zip(s_hist, y_hist, rho_hist, reversed(alphas)):
            b = rho_i * _dot(y_i, r)
            r = [rj + (a - b) * sj for rj, sj in zip(r, s_i)]
        p = [-rj for rj in r]

        slope = _dot(g, p)
        if slope >= 0:
            # not a descent direction (history stale) -> steepest descent, reset memory
            p = [-gi for gi in g]
            slope = _dot(g, p)
            s_hist.clear()
            y_hist.clear()
            rho_hist.clear()

        # ---- strong-Wolfe line search (bracket + zoom) ----
        alpha = _wolfe_line_search(func, gfun, x, fx, g, p, slope)
        x_new = [xi + alpha * pi for xi, pi in zip(x, p)]
        f_new = func(x_new)
        g_new = gfun(x_new)

        s = [x_new[i] - x[i] for i in range(n)]
        y = [g_new[i] - g[i] for i in range(n)]
        sy = _dot(s, y)
        if sy > 1e-12:                     # curvature condition -> keep pair
            s_hist.append(s)
            y_hist.append(y)
            rho_hist.append(1.0 / sy)
            if len(s_hist) > m:
                s_hist.pop(0)
                y_hist.pop(0)
                rho_hist.pop(0)

        x, fx, g = x_new, f_new, g_new

    gnorm = math.sqrt(_dot(g, g))
    return {"x": x, "fun": fx, "n_iter": n_iter, "converged": gnorm < tol,
            "grad_norm": gnorm}


def _wolfe_line_search(func, gfun, x, f0, g0, p, slope0,
                       c1=1e-4, c2=0.9, max_ls=50):
    # strong-Wolfe: sufficient decrease (c1) + curvature (c2), via bracket-and-zoom.
    n = len(x)

    def phi(a):
        xa = [x[i] + a * p[i] for i in range(n)]
        return func(xa)

    def dphi(a):
        xa = [x[i] + a * p[i] for i in range(n)]
        return _dot(gfun(xa), p)

    a_prev = 0.0
    f_prev = f0
    a = 1.0
    for i in range(max_ls):
        fa = phi(a)
        if fa > f0 + c1 * a * slope0 or (i > 0 and fa >= f_prev):
            return _zoom(phi, dphi, a_prev, a, f0, slope0, c1, c2)
        da = dphi(a)
        if abs(da) <= -c2 * slope0:
            return a
        if da >= 0:
            return _zoom(phi, dphi, a, a_prev, f0, slope0, c1, c2)
        a_prev, f_prev = a, fa
        a *= 2.0
    return a


def _zoom(phi, dphi, a_lo, a_hi, f0, slope0, c1, c2, max_zoom=50):
    for _ in range(max_zoom):
        a = 0.5 * (a_lo + a_hi)
        fa = phi(a)
        if fa > f0 + c1 * a * slope0 or fa >= phi(a_lo):
            a_hi = a
        else:
            da = dphi(a)
            if abs(da) <= -c2 * slope0:
                return a
            if da * (a_hi - a_lo) >= 0:
                a_hi = a_lo
            a_lo = a
        if abs(a_hi - a_lo) < 1e-14:
            break
    return 0.5 * (a_lo + a_hi)
