"""Numerical differentiation: gradient, Hessian, and Jacobian.

Central finite differences for a scalar function ``f(x)`` of a vector ``x`` (and a
vector function for the Jacobian). Central differences have O(h^2) error, far more
accurate than one-sided; the default step scales with the argument magnitude.
Used for Greeks, calibration gradients, and Newton steps where no analytic
derivative is available. Pure standard library.
"""


def _step(xi, rel):
    return rel * max(abs(xi), 1.0)


def gradient(f, x, rel_step=1e-6):
    """Central-difference gradient of a scalar ``f`` at the vector ``x``.

    Returns a list of ``df/dx_i``. Each partial uses ``(f(x+h e_i) - f(x-h e_i))
    / (2 h)`` with a per-coordinate step scaled by the argument magnitude.
    """
    n = len(x)
    g = [0.0] * n
    for i in range(n):
        h = _step(x[i], rel_step)
        xp = list(x); xp[i] += h
        xm = list(x); xm[i] -= h
        g[i] = (f(xp) - f(xm)) / (2.0 * h)
    return g


def hessian(f, x, rel_step=1e-4):
    """Central-difference Hessian of a scalar ``f`` at ``x`` (symmetric).

    Diagonal terms use the second-difference stencil; off-diagonal terms the
    four-point cross stencil. A larger default step than the gradient keeps the
    second differences well-conditioned. Returns an ``n x n`` list of lists.
    """
    n = len(x)
    H = [[0.0] * n for _ in range(n)]
    f0 = f(list(x))
    hs = [_step(x[i], rel_step) for i in range(n)]
    for i in range(n):
        xp = list(x); xp[i] += hs[i]
        xm = list(x); xm[i] -= hs[i]
        H[i][i] = (f(xp) - 2.0 * f0 + f(xm)) / (hs[i] * hs[i])
    for i in range(n):
        for j in range(i + 1, n):
            xpp = list(x); xpp[i] += hs[i]; xpp[j] += hs[j]
            xpm = list(x); xpm[i] += hs[i]; xpm[j] -= hs[j]
            xmp = list(x); xmp[i] -= hs[i]; xmp[j] += hs[j]
            xmm = list(x); xmm[i] -= hs[i]; xmm[j] -= hs[j]
            val = (f(xpp) - f(xpm) - f(xmp) + f(xmm)) / (4.0 * hs[i] * hs[j])
            H[i][j] = H[j][i] = val
    return H


def jacobian(f, x, rel_step=1e-6):
    """Central-difference Jacobian of a vector function ``f: R^n -> R^m``.

    Returns an ``m x n`` list of lists ``df_k/dx_i``.
    """
    n = len(x)
    f0 = f(list(x))
    m = len(f0)
    J = [[0.0] * n for _ in range(m)]
    for i in range(n):
        h = _step(x[i], rel_step)
        xp = list(x); xp[i] += h
        xm = list(x); xm[i] -= h
        fp = f(xp)
        fm = f(xm)
        for k in range(m):
            J[k][i] = (fp[k] - fm[k]) / (2.0 * h)
    return J
