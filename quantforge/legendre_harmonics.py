"""Legendre polynomials, associated Legendre functions, and real spherical harmonics.

The Legendre polynomials ``P_l(x)`` are the orthogonal polynomials on ``[-1, 1]`` behind
Gauss-Legendre quadrature, multipole expansions, and the solution of Laplace's equation. Their
associated functions ``P_l^m(x)`` and the spherical harmonics ``Y_l^m(theta, phi)`` form the
angular basis for the hydrogen atom, geopotential/geomagnetic models, and lighting in graphics.

This module evaluates ``P_l`` by the three-term recurrence, ``P_l^m`` by the stable
upward recurrence seeded from ``P_m^m`` (Numerical Recipes 6.7), and the *real* (tesseral)
spherical harmonics with full Condon-Shortley normalization so that they are orthonormal over
the sphere. Pure standard library.
"""

import math


def legendre_p(l, x):
    """Legendre polynomial ``P_l(x)`` via the three-term recurrence.

    ``P_0 = 1``, ``P_1 = x``, ``(l+1) P_{l+1} = (2l+1) x P_l - l P_{l-1}``. Valid for any real
    ``x`` (orthogonal on ``[-1, 1]``). ``P_l(1) = 1`` and ``P_l(-1) = (-1)^l``.
    """
    if l < 0:
        raise ValueError("legendre_p requires l >= 0")
    if l == 0:
        return 1.0
    p_prev = 1.0
    p_cur = x
    for k in range(1, l):
        p_next = ((2 * k + 1) * x * p_cur - k * p_prev) / (k + 1)
        p_prev, p_cur = p_cur, p_next
    return p_cur


def assoc_legendre(l, m, x):
    """Associated Legendre function ``P_l^m(x)`` for ``0 <= m <= l`` and ``|x| <= 1``.

    Uses the Condon-Shortley phase ``(-1)^m``. Seeds ``P_m^m = (-1)^m (2m-1)!! (1-x^2)^{m/2}``
    then climbs in ``l`` with ``(l-m) P_l^m = x (2l-1) P_{l-1}^m - (l+m-1) P_{l-2}^m``.
    ``assoc_legendre(l, 0, x) == legendre_p(l, x)``.
    """
    if m < 0 or m > l:
        raise ValueError("assoc_legendre requires 0 <= m <= l")
    if abs(x) > 1.0:
        raise ValueError("assoc_legendre requires |x| <= 1")
    # P_m^m
    pmm = 1.0
    if m > 0:
        somx2 = math.sqrt((1.0 - x) * (1.0 + x))
        fact = 1.0
        for _ in range(m):
            pmm *= -fact * somx2   # (-1)^m (2k-1)!! factor
            fact += 2.0
    if l == m:
        return pmm
    # P_{m+1}^m
    pmmp1 = x * (2 * m + 1) * pmm
    if l == m + 1:
        return pmmp1
    # climb to P_l^m
    pll = 0.0
    for ll in range(m + 2, l + 1):
        pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m)
        pmm = pmmp1
        pmmp1 = pll
    return pll


def _norm(l, m):
    # sqrt((2l+1)/(4pi) * (l-m)!/(l+m)!)
    num = (2 * l + 1)
    ratio = 1.0
    for k in range(l - m + 1, l + m + 1):
        ratio *= k          # (l+m)!/(l-m)! = product
    return math.sqrt(num / (4.0 * math.pi * ratio))


def spherical_harmonic_real(l, m, theta, phi):
    """Real spherical harmonic ``Y_l^m(theta, phi)`` (orthonormal, ``-l <= m <= l``).

    ``theta`` is the polar (colatitude) angle in ``[0, pi]``, ``phi`` the azimuth. The real
    (tesseral) convention:
    ``m > 0``: ``sqrt(2) N_l^m P_l^m(cos theta) cos(m phi)``;
    ``m = 0``: ``N_l^0 P_l^0(cos theta)``;
    ``m < 0``: ``sqrt(2) N_l^|m| P_l^|m|(cos theta) sin(|m| phi)``,
    with ``N_l^m = sqrt((2l+1)/(4 pi) (l-m)!/(l+m)!)``. These are orthonormal over the sphere.
    """
    if l < 0 or abs(m) > l:
        raise ValueError("spherical_harmonic_real requires 0 <= |m| <= l")
    am = abs(m)
    p = assoc_legendre(l, am, math.cos(theta))
    n = _norm(l, am)
    if m == 0:
        return n * p
    if m > 0:
        return math.sqrt(2.0) * n * p * math.cos(m * phi)
    return math.sqrt(2.0) * n * p * math.sin(am * phi)
