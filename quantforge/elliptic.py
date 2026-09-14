"""Complete elliptic integrals K and E, and the arithmetic-geometric mean.

The complete elliptic integrals appear in pendulum periods, ellipse arc length, and
electromagnetics. ``K(m) = integral_0^{pi/2} dtheta / sqrt(1 - m sin^2 theta)`` (first kind)
and ``E(m) = integral_0^{pi/2} sqrt(1 - m sin^2 theta) dtheta`` (second kind), with parameter
``m = k^2`` in ``[0, 1)`` (and ``m <= 1`` for E). Both are evaluated from the
arithmetic-geometric mean, which converges quadratically: ``K = pi / (2 AGM(1, sqrt(1-m)))``.
Pure standard library.
"""

import math


def agm(a, b, tol=1e-15, max_iter=100):
    """Arithmetic-geometric mean of ``a, b >= 0``: the common limit of the AM/GM iteration.

    Each step replaces ``(a, b)`` with ``((a+b)/2, sqrt(a b))``; the two converge
    quadratically. ``agm(a, b) == agm(b, a)`` and ``agm(a, a) == a``.
    """
    if a < 0 or b < 0:
        raise ValueError("agm requires non-negative arguments")
    for _ in range(max_iter):
        a2 = 0.5 * (a + b)
        b2 = math.sqrt(a * b)
        if abs(a2 - b2) <= tol * a2 if a2 else abs(a2 - b2) <= tol:
            return a2
        a, b = a2, b2
    return 0.5 * (a + b)


def elliptic_k(m):
    """Complete elliptic integral of the first kind ``K(m)``, ``0 <= m < 1``.

    ``K(0) = pi/2`` and ``K(m) -> inf`` as ``m -> 1``. Computed as
    ``pi / (2 AGM(1, sqrt(1 - m)))``.
    """
    if not 0 <= m < 1:
        raise ValueError("elliptic_k requires 0 <= m < 1")
    return math.pi / (2 * agm(1.0, math.sqrt(1 - m)))


def elliptic_e(m):
    """Complete elliptic integral of the second kind ``E(m)``, ``0 <= m <= 1``.

    ``E(0) = pi/2`` and ``E(1) = 1``. Uses the AGM descent, accumulating the geometric-step
    corrections that give ``E`` from the same iteration.
    """
    if not 0 <= m <= 1:
        raise ValueError("elliptic_e requires 0 <= m <= 1")
    if m == 1.0:
        return 1.0
    a = 1.0
    b = math.sqrt(1 - m)
    c = math.sqrt(m)              # c_0 = sqrt(a0^2 - b0^2)
    s = 0.5 * c * c               # running sum: sum_{n>=0} 2^(n-1) c_n^2, n=0 term
    powr = 1.0                    # weight 2^(n-1) for the next c_n; c_1 gets 2^0 = 1
    for _ in range(100):
        a2 = 0.5 * (a + b)
        b2 = math.sqrt(a * b)
        c = 0.5 * (a - b)         # c_{n+1}
        s += powr * c * c
        powr *= 2.0
        if abs(a - b) <= 1e-15 * a:
            a = a2
            break
        a, b = a2, b2
    k = math.pi / (2 * a)
    return k * (1 - s)
