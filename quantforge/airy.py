"""Airy functions Ai(x) and Bi(x).

The Airy functions solve ``y'' = x y`` -- the turning-point equation of quantum mechanics
and optics (caustics, the rainbow). ``Ai`` decays for ``x > 0``; ``Bi`` grows. Both are
entire, so they are built from the two power-series solutions ``f`` and ``g`` of the ODE:
``Ai = c1 f - c2 g`` and ``Bi = sqrt(3)(c1 f + c2 g)`` with the standard constants. The
series converges everywhere; for large positive ``x`` the terms are well-behaved because the
combination decays. Note: for large *negative* ``x`` (roughly ``x < -15``) the two series
terms nearly cancel, so accuracy degrades there -- the reliable range is about
``[-15, +inf)``. Pure standard library.
"""

import math

# Ai(0) = 1 / (3^(2/3) Gamma(2/3)); Ai'(0) = -1 / (3^(1/3) Gamma(1/3))
_C1 = 0.3550280538878172   # Ai(0)
_C2 = 0.2588194037928068   # -Ai'(0)


def _fg(x, nterms=200):
    # f(x) = sum x^{3k} / (3k)! * prod, g(x) = sum x^{3k+1}/(3k+1)! * prod
    # f: 1 + x^3/6 + ... coefficients a_{k}: a0=1, a_{k}= a_{k-1} * 1 / ((3k-1)(3k))? use recurrence
    # Standard series:
    # f = sum_{k>=0} 3^k (1/3)_k / (3k)! x^{3k}
    # g = sum_{k>=0} 3^k (2/3)_k / (3k+1)! x^{3k+1}
    f = 0.0
    g = 0.0
    # f term
    term_f = 1.0
    for k in range(nterms):
        f += term_f
        # next: multiply by x^3 * (3k+1) / ((3k+1)(3k+2)(3k+3)) = x^3/((3k+2)(3k+3))
        term_f *= x ** 3 / ((3 * k + 2) * (3 * k + 3))
        if abs(term_f) < 1e-20 * (abs(f) + 1e-30):
            break
    term_g = x
    for k in range(nterms):
        g += term_g
        term_g *= x ** 3 / ((3 * k + 3) * (3 * k + 4))
        if abs(term_g) < 1e-20 * (abs(g) + 1e-30):
            break
    return f, g


def airy_ai(x):
    """Airy function of the first kind ``Ai(x)``, solving ``y'' = x y`` with ``Ai -> 0`` as ``x -> inf``."""
    f, g = _fg(x)
    return _C1 * f - _C2 * g


def airy_bi(x):
    """Airy function of the second kind ``Bi(x)`` (the growing solution)."""
    f, g = _fg(x)
    return math.sqrt(3.0) * (_C1 * f + _C2 * g)
