"""Polylogarithm and dilogarithm functions.

The polylogarithm ``Li_s(z) = sum_{k>=1} z^k / k^s`` interpolates a whole family of classical
functions: ``Li_1(z) = -ln(1 - z)``, ``Li_2`` is the dilogarithm (Spence's function), and at
``z = 1`` every ``Li_s`` collapses to the Riemann zeta ``zeta(s)``. It shows up in Fermi-Dirac
integrals, one-loop Feynman diagrams, and entropy-of-entanglement sums.

:func:`polylog` evaluates ``Li_s(z)`` by direct summation for real ``s`` and real ``|z| < 1``
(with the exact endpoints ``z = 1 -> zeta(s)`` and ``z = -1 -> -eta(s)``). :func:`dilog`
handles the dilogarithm ``Li_2(x)`` on the entire real half-line ``x <= 1`` using the standard
reflection/inversion identities to fold every argument into the fast-converging window
``[-1, 1/2]``. Pure standard library.
"""

import math

from .zeta import riemann_zeta, dirichlet_eta


def polylog(s, z, tol=1e-15, max_terms=200000):
    """Polylogarithm ``Li_s(z) = sum_{k>=1} z^k / k^s`` for real ``s`` and real ``|z| <= 1``.

    Uses direct summation, which converges geometrically for ``|z| < 1``. The endpoint
    ``z = 1`` returns ``zeta(s)`` (requires ``s > 1``); ``z = -1`` returns ``-eta(s)``.
    Convergence slows as ``z -> 1`` with small ``s``; ``max_terms`` bounds the work.
    """
    if z > 1.0 or z < -1.0:
        raise ValueError("polylog supports real |z| <= 1")
    if z == 1.0:
        return riemann_zeta(s)          # raises for s <= 1
    if z == -1.0:
        return -dirichlet_eta(s)        # Li_s(-1) = -eta(s), s > 0
    if z == 0.0:
        return 0.0
    total = 0.0
    zk = z
    for k in range(1, max_terms + 1):
        term = zk / (k ** s)
        total += term
        if abs(term) < tol * (abs(total) + tol):
            break
        zk *= z
    return total


def dilog(x):
    """Dilogarithm (Spence's function) ``Li_2(x) = sum_{k>=1} x^k / k^2`` for real ``x <= 1``.

    Reflection and inversion identities fold ``x`` into ``[-1, 1/2]`` where the series
    converges quickly. ``Li_2(1) = pi^2/6``, ``Li_2(-1) = -pi^2/12``, and
    ``Li_2(1/2) = pi^2/12 - (ln 2)^2/2``. Raises for ``x > 1`` (there ``Li_2`` is complex).
    """
    if x > 1.0:
        raise ValueError("dilog is real only for x <= 1")
    pi2_6 = math.pi ** 2 / 6.0
    if x == 1.0:
        return pi2_6
    if x == 0.0:
        return 0.0
    # Fold large-magnitude arguments inward.
    if x < -1.0:
        # inversion: Li_2(x) = -Li_2(1/x) - pi^2/6 - (1/2) ln(-x)^2
        return -_dilog_series(1.0 / x) - pi2_6 - 0.5 * math.log(-x) ** 2
    if x > 0.5:
        # reflection: Li_2(x) = pi^2/6 - ln(x) ln(1-x) - Li_2(1-x)
        return pi2_6 - math.log(x) * math.log(1.0 - x) - _dilog_series(1.0 - x)
    return _dilog_series(x)


def _dilog_series(x, tol=1e-17, max_terms=100000):
    # direct series, |x| <= 1/2 after folding -> fast geometric convergence
    total = 0.0
    xk = x
    for k in range(1, max_terms + 1):
        term = xk / (k * k)
        total += term
        if abs(term) < tol * (abs(total) + tol):
            break
        xk *= x
    return total
