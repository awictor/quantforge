"""Carlson symmetric elliptic integrals and the incomplete elliptic integrals.

Carlson's symmetric forms ``R_F, R_C, R_D, R_J`` are the modern basis for elliptic integrals:
they are numerically stable, symmetric in their arguments, and evaluate every incomplete
(Legendre) integral through simple algebraic combinations. Each is computed by the
duplication-theorem iteration (Carlson 1979, Numerical Recipes 6.11), which drives the
arguments together quadratically and finishes with a short Taylor series.

From them come the incomplete elliptic integrals of the first, second, and third kind:
``F(phi | m) = sin(phi) R_F(cos^2, 1 - m sin^2, 1)`` and companions. Parameter convention is
``m = k^2`` (matching :mod:`quantforge.elliptic`). Pure standard library.
"""

import math

_TINY = 1e-300


def carlson_rf(x, y, z, tol=1e-12):
    """Carlson's ``R_F(x, y, z) = 1/2 integral_0^inf dt / sqrt((t+x)(t+y)(t+z))``.

    Symmetric and homogeneous of degree ``-1/2``. At most one argument may be zero; all must be
    non-negative. ``R_F(x, x, x) = 1/sqrt(x)``.
    """
    if x < 0 or y < 0 or z < 0:
        raise ValueError("carlson_rf requires non-negative arguments")
    if (x + y) < _TINY or (x + z) < _TINY or (y + z) < _TINY:
        raise ValueError("carlson_rf: at most one argument may be zero")
    for _ in range(100):
        sx = math.sqrt(x)
        sy = math.sqrt(y)
        sz = math.sqrt(z)
        lam = sx * sy + sy * sz + sz * sx
        x = 0.25 * (x + lam)
        y = 0.25 * (y + lam)
        z = 0.25 * (z + lam)
        mu = (x + y + z) / 3.0
        dx = 1.0 - x / mu
        dy = 1.0 - y / mu
        dz = 1.0 - z / mu
        if max(abs(dx), abs(dy), abs(dz)) < tol:
            break
    e2 = dx * dy - dz * dz
    e3 = dx * dy * dz
    # Taylor series (NR 6.11.10)
    return (1.0 + (e2 / 24.0 - 0.1 - 3.0 * e3 / 44.0) * e2
            + e3 / 14.0) / math.sqrt(mu)


def carlson_rc(x, y, tol=1e-12):
    """Carlson's degenerate ``R_C(x, y) = R_F(x, y, y)``.

    ``R_C(x, x) = 1/sqrt(x)``; equals ``arctan``/``arctanh``-type elementary functions.
    """
    return carlson_rf(x, y, y, tol)


def carlson_rd(x, y, z, tol=1e-12):
    """Carlson's ``R_D(x, y, z) = R_J(x, y, z, z)``, symmetric in ``x, y`` only.

    ``= 3/2 integral_0^inf dt / [(t+z) sqrt((t+x)(t+y)(t+z))]``. ``z`` must be positive; ``x, y``
    non-negative with at most one zero. Homogeneous of degree ``-3/2``.
    """
    if x < 0 or y < 0:
        raise ValueError("carlson_rd requires non-negative x, y")
    if z <= 0:
        raise ValueError("carlson_rd requires positive z")
    total = 0.0
    fac = 1.0
    for _ in range(100):
        sx = math.sqrt(x)
        sy = math.sqrt(y)
        sz = math.sqrt(z)
        lam = sx * sy + sy * sz + sz * sx
        total += fac / (sz * (z + lam))
        fac *= 0.25
        x = 0.25 * (x + lam)
        y = 0.25 * (y + lam)
        z = 0.25 * (z + lam)
        mu = (x + y + 3.0 * z) / 5.0
        dx = 1.0 - x / mu
        dy = 1.0 - y / mu
        dz = 1.0 - z / mu
        if max(abs(dx), abs(dy), abs(dz)) < tol:
            break
    ea = dx * dy
    eb = dz * dz
    ec = ea - eb
    ed = ea - 6.0 * eb
    ee = ed + ec + ec
    # NR 6.11.13
    s = (ed * (-3.0 / 14.0 + 0.25 * ed - 1.5 * dz * ee / 11.0)
         + dz * (ee / 6.0 + dz * (-ec * 9.0 / 22.0 + dz * ea * 3.0 / 26.0)))
    return 3.0 * total + fac * (1.0 + s) / (mu * math.sqrt(mu))


def carlson_rj(x, y, z, p, tol=1e-12):
    """Carlson's ``R_J(x, y, z, p) = 3/2 integral_0^inf dt / [(t+p) sqrt((t+x)(t+y)(t+z))]``.

    Symmetric in ``x, y, z``; all non-negative with at most one zero, and ``p != 0``. This
    implementation covers ``p > 0``.
    """
    if x < 0 or y < 0 or z < 0:
        raise ValueError("carlson_rj requires non-negative x, y, z")
    if p <= 0:
        raise ValueError("carlson_rj implemented for p > 0")
    total = 0.0
    fac = 1.0
    x0, y0, z0, p0 = x, y, z, p
    for _ in range(100):
        sx = math.sqrt(x)
        sy = math.sqrt(y)
        sz = math.sqrt(z)
        lam = sx * sy + sy * sz + sz * sx
        alpha = (p * (sx + sy + sz) + sx * sy * sz) ** 2
        beta = p * (p + lam) ** 2
        total += fac * carlson_rc(alpha, beta, tol)
        fac *= 0.25
        x = 0.25 * (x + lam)
        y = 0.25 * (y + lam)
        z = 0.25 * (z + lam)
        p = 0.25 * (p + lam)
        mu = (x + y + z + 2.0 * p) / 5.0
        dx = 1.0 - x / mu
        dy = 1.0 - y / mu
        dz = 1.0 - z / mu
        dp = 1.0 - p / mu
        if max(abs(dx), abs(dy), abs(dz), abs(dp)) < tol:
            break
    c1 = 3.0 / 14.0
    c2 = 1.0 / 3.0
    c3 = 3.0 / 22.0
    c4 = 3.0 / 26.0
    ea = dx * (dy + dz) + dy * dz
    eb = dx * dy * dz
    ec = dp * dp
    ed = ea - 3.0 * ec
    ee = eb + 2.0 * dp * (ea - ec)
    # NR 6.11.20
    s = (1.0 + ed * (-c1 + 0.25 * c3 * ed - 1.5 * c4 * ee)
         + eb * (c2 + dp * (-c3 - c3 + dp * c4))
         + dp * ea * (c2 - dp * c3) - c2 * dp * ec)
    return 3.0 * total + fac * s / (mu * math.sqrt(mu))


def elliptic_f(phi, m, tol=1e-12):
    """Incomplete elliptic integral of the first kind ``F(phi | m)`` via ``R_F``.

    ``F(phi | m) = integral_0^phi dtheta / sqrt(1 - m sin^2 theta) = sin(phi) R_F(c, 1 - m s^2, 1)``
    with ``s = sin phi``, ``c = cos^2 phi``. ``m = k^2``; reduces to ``F(pi/2, m) = K(m)``.
    """
    s = math.sin(phi)
    c = math.cos(phi)
    return s * carlson_rf(c * c, 1.0 - m * s * s, 1.0, tol)


def elliptic_e_incomplete(phi, m, tol=1e-12):
    """Incomplete elliptic integral of the second kind ``E(phi | m)`` via ``R_F`` and ``R_D``.

    ``E(phi | m) = sin(phi) R_F(c, 1-m s^2, 1) - (m/3) sin^3(phi) R_D(c, 1-m s^2, 1)``.
    Reduces to the complete ``E(m)`` at ``phi = pi/2``.
    """
    s = math.sin(phi)
    c = math.cos(phi)
    cc = c * c
    q = 1.0 - m * s * s
    return (s * carlson_rf(cc, q, 1.0, tol)
            - (m / 3.0) * s ** 3 * carlson_rd(cc, q, 1.0, tol))


def elliptic_pi(n, phi, m, tol=1e-12):
    """Incomplete elliptic integral of the third kind ``Pi(n; phi | m)`` via ``R_F`` and ``R_J``.

    ``Pi(n; phi | m) = integral_0^phi dtheta / [(1 - n sin^2 theta) sqrt(1 - m sin^2 theta)]``,
    ``= sin(phi) R_F(c, 1-m s^2, 1) + (n/3) sin^3(phi) R_J(c, 1-m s^2, 1, 1 - n s^2)``.
    Implemented for ``n < 1`` (so the ``R_J`` fourth argument stays positive).
    """
    s = math.sin(phi)
    c = math.cos(phi)
    cc = c * c
    q = 1.0 - m * s * s
    p = 1.0 - n * s * s
    return (s * carlson_rf(cc, q, 1.0, tol)
            + (n / 3.0) * s ** 3 * carlson_rj(cc, q, 1.0, p, tol))
