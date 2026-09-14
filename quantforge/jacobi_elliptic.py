"""Jacobi elliptic functions sn, cn, dn and the amplitude am.

The Jacobi elliptic functions invert the incomplete elliptic integral of the first kind: if
``u = integral_0^phi dtheta / sqrt(1 - m sin^2 theta)`` then ``am(u, m) = phi`` and
``sn(u, m) = sin(phi)``, ``cn(u, m) = cos(phi)``, ``dn(u, m) = sqrt(1 - m sin^2 phi)``. They are
the doubly-periodic generalization of sine/cosine (recovered at ``m = 0``) and describe the
exact pendulum, solitons, and the motion of a rigid body. This uses the descending-Landen /
arithmetic-geometric-mean algorithm (Abramowitz & Stegun 16.4), which converges quadratically.
Pure standard library.
"""

import math

from .elliptic import elliptic_k


def jacobi_am(u, m, tol=1e-15, max_iter=64):
    """Jacobi amplitude ``am(u, m)`` -- the angle ``phi`` with ``u = F(phi | m)``.

    AGM descent: build the sequence ``a_0 = 1, b_0 = sqrt(1-m), c_0 = sqrt(m)`` down to
    ``c_n ~ 0``, then climb back, halving the accumulated angle. ``m`` is the parameter
    ``k^2`` in ``[0, 1]``. Reduces to ``am(u, 0) = u``.
    """
    if not 0.0 <= m <= 1.0:
        raise ValueError("jacobi_am requires 0 <= m <= 1")
    if m == 0.0:
        return u
    if m == 1.0:
        # degenerate: am(u, 1) = 2 arctan(e^u) - pi/2 (gudermannian)
        return 2.0 * math.atan(math.exp(u)) - math.pi / 2.0
    a = 1.0
    b = math.sqrt(1.0 - m)
    c = math.sqrt(m)
    a_list = [a]
    c_list = [c]
    n = 0
    while abs(c) > tol * abs(a) and n < max_iter:
        a2 = 0.5 * (a + b)
        b2 = math.sqrt(a * b)
        c = 0.5 * (a - b)
        a, b = a2, b2
        a_list.append(a)
        c_list.append(c)
        n += 1
    # phi_n = 2^n * a_n * u
    phi = (2 ** n) * a * u
    # descend the angle chain
    for i in range(n, 0, -1):
        phi = 0.5 * (phi + math.asin(c_list[i] / a_list[i] * math.sin(phi)))
    return phi


def jacobi_sn(u, m):
    """Jacobi elliptic ``sn(u, m) = sin(am(u, m))``. Reduces to ``sin u`` at ``m = 0``."""
    return math.sin(jacobi_am(u, m))


def jacobi_cn(u, m):
    """Jacobi elliptic ``cn(u, m) = cos(am(u, m))``. Reduces to ``cos u`` at ``m = 0``."""
    return math.cos(jacobi_am(u, m))


def jacobi_dn(u, m):
    """Jacobi elliptic ``dn(u, m) = sqrt(1 - m sin^2 am(u, m))``. Reduces to ``1`` at ``m = 0``."""
    s = math.sin(jacobi_am(u, m))
    return math.sqrt(1.0 - m * s * s)
