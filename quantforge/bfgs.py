"""BFGS quasi-Newton minimization with a backtracking line search.

For a smooth objective, gradient methods converge far faster than the derivative-free
simplex or population searches. Newton's method needs the Hessian; BFGS (Broyden-
Fletcher-Goldfarb-Shanno) instead builds an approximation to the *inverse* Hessian from
successive gradients, giving near-Newton (superlinear) convergence using only gradients.
The gradient is taken numerically by central differences, so only the objective is
required. An Armijo backtracking line search guarantees descent. Pure standard library.
"""

import math


def _grad(f, x, h=1e-6):
    n = len(x)
    g = [0.0] * n
    for i in range(n):
        step = h * (abs(x[i]) + 1.0)
        xp = list(x); xp[i] += step
        xm = list(x); xm[i] -= step
        g[i] = (f(xp) - f(xm)) / (2.0 * step)
    return g


def bfgs(func, x0, tol=1e-8, max_iter=500):
    """Minimize ``func`` from ``x0`` by BFGS with a backtracking line search.

    ``func`` takes a length-``n`` list and returns a scalar; the gradient is computed by
    central differences. Returns a dict with ``x`` (minimizer), ``fun`` (its value),
    ``n_iter``, ``converged`` (gradient norm below ``tol``) and ``grad_norm``. Best for
    smooth objectives; use a global method first if the landscape is multimodal.
    """
    n = len(x0)
    if n == 0:
        raise ValueError("need at least one dimension")
    x = [float(v) for v in x0]
    g = _grad(func, x)
    # Inverse-Hessian approximation, initialized to the identity.
    H = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    fx = func(x)

    n_iter = 0
    for it in range(max_iter):
        n_iter = it + 1
        gnorm = math.sqrt(sum(gi * gi for gi in g))
        if gnorm < tol:
            return {"x": x, "fun": fx, "n_iter": n_iter, "converged": True,
                    "grad_norm": gnorm}
        # Search direction p = -H g.
        p = [-sum(H[i][j] * g[j] for j in range(n)) for i in range(n)]
        # Armijo backtracking line search.
        slope = sum(g[i] * p[i] for i in range(n))
        if slope >= 0:
            # H lost positive-definiteness; reset to steepest descent.
            p = [-gi for gi in g]
            slope = sum(g[i] * p[i] for i in range(n))
            H = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
        alpha = 1.0
        c1 = 1e-4
        while True:
            x_new = [x[i] + alpha * p[i] for i in range(n)]
            f_new = func(x_new)
            if f_new <= fx + c1 * alpha * slope or alpha < 1e-14:
                break
            alpha *= 0.5
        g_new = _grad(func, x_new)
        # BFGS update of H using s = x_new - x, y = g_new - g.
        s = [x_new[i] - x[i] for i in range(n)]
        y = [g_new[i] - g[i] for i in range(n)]
        sy = sum(s[i] * y[i] for i in range(n))
        if sy > 1e-12:
            rho = 1.0 / sy
            # H <- (I - rho s y') H (I - rho y s') + rho s s'.
            Hy = [sum(H[i][j] * y[j] for j in range(n)) for i in range(n)]
            yHy = sum(y[i] * Hy[i] for i in range(n))
            newH = [[H[i][j]
                     + rho * rho * yHy * s[i] * s[j]
                     - rho * (s[i] * Hy[j] + Hy[i] * s[j])
                     + rho * s[i] * s[j]
                     for j in range(n)] for i in range(n)]
            H = newH
        x, g, fx = x_new, g_new, f_new

    gnorm = math.sqrt(sum(gi * gi for gi in g))
    return {"x": x, "fun": fx, "n_iter": n_iter, "converged": gnorm < tol,
            "grad_norm": gnorm}
