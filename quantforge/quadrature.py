"""Numerical integration: trapezoid, Simpson, Gauss-Legendre, adaptive Simpson,
Romberg, and tanh-sinh.

General-purpose definite-integral routines for the many expected-value and
density integrals across the library. Trapezoid and composite Simpson on a fixed
grid, fixed-order Gauss-Legendre (exact for polynomials up to degree 2n-1), an
error-controlled adaptive Simpson, Romberg (Richardson extrapolation on the
trapezoid rule) for smooth integrands, and tanh-sinh (double-exponential) for
integrable endpoint singularities. Pure standard library.
"""

import math


def trapezoid(f, a, b, n=1000):
    """Composite trapezoid rule with ``n`` sub-intervals over ``[a, b]``."""
    if n < 1:
        raise ValueError("n must be positive")
    h = (b - a) / n
    total = 0.5 * (f(a) + f(b))
    for i in range(1, n):
        total += f(a + i * h)
    return total * h


def simpson(f, a, b, n=1000):
    """Composite Simpson's rule (``n`` even) -- exact for cubics.

    Rounds ``n`` up to the next even number. Fourth-order accurate.
    """
    if n < 2:
        raise ValueError("n must be >= 2")
    if n % 2:
        n += 1
    h = (b - a) / n
    total = f(a) + f(b)
    for i in range(1, n):
        total += (4.0 if i % 2 else 2.0) * f(a + i * h)
    return total * h / 3.0


# Gauss-Legendre nodes/weights on [-1, 1] for a few common orders.
_GL = {
    2: ([-0.5773502691896257, 0.5773502691896257], [1.0, 1.0]),
    3: ([-0.7745966692414834, 0.0, 0.7745966692414834],
        [0.5555555555555556, 0.8888888888888888, 0.5555555555555556]),
    4: ([-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
         0.8611363115940526],
        [0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
         0.3478548451374538]),
    5: ([-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831,
         0.9061798459386640],
        [0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
         0.4786286704993665, 0.2369268850561891]),
}


def gauss_legendre(f, a, b, n=5):
    """Fixed-order Gauss-Legendre quadrature (``n`` in {2,3,4,5}).

    Maps the reference nodes to ``[a, b]``. Exact for polynomials up to degree
    ``2n - 1`` -- very accurate for smooth integrands with few evaluations.
    """
    if n not in _GL:
        raise ValueError("n must be one of 2, 3, 4, 5")
    nodes, weights = _GL[n]
    half = 0.5 * (b - a)
    mid = 0.5 * (a + b)
    return half * sum(weights[i] * f(mid + half * nodes[i]) for i in range(n))


def _adaptive(f, a, b, fa, fb, fm, whole, tol, depth):
    m = 0.5 * (a + b)
    lm = 0.5 * (a + m)
    rm = 0.5 * (m + b)
    flm = f(lm)
    frm = f(rm)
    left = (m - a) / 6.0 * (fa + 4.0 * flm + fm)
    right = (b - m) / 6.0 * (fm + 4.0 * frm + fb)
    if depth <= 0 or abs(left + right - whole) <= 15.0 * tol:
        return left + right + (left + right - whole) / 15.0
    return (_adaptive(f, a, m, fa, fm, flm, left, tol / 2.0, depth - 1)
            + _adaptive(f, m, b, fm, fb, frm, right, tol / 2.0, depth - 1))


def adaptive_simpson(f, a, b, tol=1e-10, max_depth=50):
    """Adaptive Simpson quadrature with error control to ``tol``.

    Recursively bisects where the Simpson estimate has not converged, so it
    concentrates work on the hard parts of the integrand. Returns the integral.
    """
    m = 0.5 * (a + b)
    fa, fb, fm = f(a), f(b), f(m)
    whole = (b - a) / 6.0 * (fa + 4.0 * fm + fb)
    return _adaptive(f, a, b, fa, fb, fm, whole, tol, max_depth)


def romberg(f, a, b, max_order=10, tol=1e-12):
    """Romberg integration: Richardson extrapolation on the trapezoid rule.

    Builds the Romberg tableau, refining the composite trapezoid estimate by
    successive interval halvings and extrapolating away the Euler-Maclaurin error
    terms. Row ``T[k][0]`` is the ``2^k``-panel trapezoid; each further column
    cancels the next even power of the step,

        T[k][j] = (4^j T[k][j-1] - T[k-1][j-1]) / (4^j - 1),

    so ``T[k][k]`` converges as ``O(h^{2k+2})`` for a smooth integrand. Stops early
    when two successive diagonal estimates agree to ``tol``. Ideal for smooth
    integrands where it reaches machine precision in a handful of halvings; for
    endpoint singularities use :func:`tanh_sinh` instead.
    """
    if max_order < 1:
        raise ValueError("max_order must be at least 1")
    table = [[0.0] * (max_order + 1) for _ in range(max_order + 1)]
    h = b - a
    table[0][0] = 0.5 * h * (f(a) + f(b))
    prev_diag = table[0][0]
    for k in range(1, max_order + 1):
        h *= 0.5
        # New midpoints added at this refinement level.
        total = 0.0
        for i in range(1, 2 ** k, 2):
            total += f(a + i * h)
        table[k][0] = 0.5 * table[k - 1][0] + h * total
        pow4 = 1.0
        for j in range(1, k + 1):
            pow4 *= 4.0
            table[k][j] = (pow4 * table[k][j - 1] - table[k - 1][j - 1]) / (pow4 - 1.0)
        if abs(table[k][k] - prev_diag) <= tol:
            return table[k][k]
        prev_diag = table[k][k]
    return table[max_order][max_order]


def tanh_sinh(f, a, b, levels=6, h0=1.0):
    """Tanh-sinh (double-exponential) quadrature over ``[a, b]``.

    Substitutes ``x = (a+b)/2 + (b-a)/2 * tanh((pi/2) sinh(t))`` and integrates the
    transformed integrand in ``t`` on a uniform grid. The change of variables makes
    the abscissae cluster double-exponentially toward the endpoints and the weights
    decay super-fast, so the rule converges even when ``f`` has integrable
    endpoint singularities (``1/sqrt(x)``, ``ln x``, ...) where Gauss-Legendre and
    Simpson struggle. ``levels`` successive grid halvings refine the step from
    ``h0``; the number of function evaluations is ``O(2^levels / h0)``.

    Nodes are computed off the endpoints, so ``f`` is never evaluated exactly at
    ``a`` or ``b`` -- an integrable singularity at either end is fine.
    """
    if levels < 0:
        raise ValueError("levels must be non-negative")
    half = 0.5 * (b - a)
    mid = 0.5 * (a + b)
    half_pi = 0.5 * math.pi

    def contrib(t):
        # Weight and abscissa of the double-exponential transform at parameter t.
        s = math.sinh(t)
        u = half_pi * s
        x = math.tanh(u)                      # in (-1, 1)
        # dx/dt = (pi/2) cosh(t) / cosh(u)^2, guarded against overflow.
        cosh_u = math.cosh(u)
        w = half_pi * math.cosh(t) / (cosh_u * cosh_u)
        return x, w

    # Finest step and total half-width in t: beyond t_max the weight underflows.
    h = h0 / (2 ** levels)
    t_max = 3.5
    n = int(t_max / h)
    total = 0.0
    for k in range(-n, n + 1):
        t = k * h
        x, w = contrib(t)
        if abs(x) >= 1.0:
            continue
        total += w * f(mid + half * x)
    return total * h * half
