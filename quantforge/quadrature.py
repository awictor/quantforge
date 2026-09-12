"""Numerical integration: trapezoid, Simpson, Gauss-Legendre, adaptive Simpson.

General-purpose definite-integral routines for the many expected-value and
density integrals across the library. Trapezoid and composite Simpson on a fixed
grid, fixed-order Gauss-Legendre (exact for polynomials up to degree 2n-1), and
an error-controlled adaptive Simpson. Pure standard library.
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
