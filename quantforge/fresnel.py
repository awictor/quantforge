"""Fresnel integrals S(x), C(x) and the Dawson function.

The Fresnel integrals ``C(x) = integral_0^x cos(pi t^2 / 2) dt`` and
``S(x) = integral_0^x sin(pi t^2 / 2) dt`` trace the Cornu spiral and describe near-field
diffraction; both approach ``1/2`` as ``x -> inf``. The Dawson function
``D(x) = e^{-x^2} integral_0^x e^{t^2} dt`` is the scaled imaginary error function, peaking
near ``x = 0.924``. Fresnel uses the convergent power series for small ``x`` and the
auxiliary-function asymptotics for large ``x``; Dawson uses its power series and continued
fraction. Pure standard library.
"""

import math


def _fresnel_series(x):
    # power series for C and S, good for |x| up to ~1.6
    s = 0.0
    c = 0.0
    t = x
    # C: sum (-1)^n (pi/2)^{2n} x^{4n+1} / ((2n)! (4n+1))
    # S: sum (-1)^n (pi/2)^{2n+1} x^{4n+3} / ((2n+1)! (4n+3))
    half_pi = math.pi / 2
    n = 0
    term_c = x
    c = x / 1.0
    # build C
    fact = 1.0
    powx = x
    acc_c = 0.0
    for n in range(0, 60):
        # C term: (-1)^n (pi/2)^{2n} x^{4n+1}/((2n)!(4n+1))
        num = ((-1) ** n) * (half_pi ** (2 * n)) * (x ** (4 * n + 1))
        den = math.factorial(2 * n) * (4 * n + 1)
        d = num / den
        acc_c += d
        if abs(d) < 1e-18 * (abs(acc_c) + 1e-30):
            break
    acc_s = 0.0
    for n in range(0, 60):
        num = ((-1) ** n) * (half_pi ** (2 * n + 1)) * (x ** (4 * n + 3))
        den = math.factorial(2 * n + 1) * (4 * n + 3)
        d = num / den
        acc_s += d
        if abs(d) < 1e-18 * (abs(acc_s) + 1e-30):
            break
    return acc_c, acc_s


def _fresnel_asymptotic(x):
    # for large x: C = 1/2 + f sin(z) - g cos(z), S = 1/2 - f cos(z) - g sin(z)
    z = math.pi / 2 * x * x
    # auxiliary functions f(x), g(x) via their asymptotic series in t = 1/(pi x^2)
    inv = 1.0 / (math.pi * x)
    t = 1.0 / (math.pi * x * x)
    f = inv * (1.0 - 3.0 * t * t + 105.0 * t ** 4)
    g = inv * t * (1.0 - 15.0 * t * t + 945.0 * t ** 4)
    sz, cz = math.sin(z), math.cos(z)
    c = 0.5 + f * sz - g * cz
    s = 0.5 - f * cz - g * sz
    return c, s


def fresnel_c(x):
    """Fresnel cosine integral ``C(x) = integral_0^x cos(pi t^2 / 2) dt``."""
    sign = 1.0 if x >= 0 else -1.0
    ax = abs(x)
    if ax < 4.0:
        c, _ = _fresnel_series(ax)
    else:
        c, _ = _fresnel_asymptotic(ax)
    return sign * c


def fresnel_s(x):
    """Fresnel sine integral ``S(x) = integral_0^x sin(pi t^2 / 2) dt``."""
    sign = 1.0 if x >= 0 else -1.0
    ax = abs(x)
    if ax < 4.0:
        _, s = _fresnel_series(ax)
    else:
        _, s = _fresnel_asymptotic(ax)
    return sign * s


def dawson(x):
    """Dawson function ``D(x) = e^{-x^2} integral_0^x e^{t^2} dt``."""
    ax = abs(x)
    if ax < 3.0:
        # power series: D(x) = sum_{n>=0} (-1)^n 2^n / (2n+1)!! x^{2n+1}
        x2 = x * x
        term = x
        s = x
        for n in range(1, 200):
            term *= -2.0 * x2 / (2 * n + 1)
            s += term
            if abs(term) < 1e-18 * (abs(s) + 1e-30):
                break
        return s
    # asymptotic series for large |x|: D ~ 1/(2x) (1 + 1/(2x^2) + 3/(2x^2)^2 + ...)
    inv = 1.0 / (2 * x)
    t = 1.0 / (2 * x * x)
    s = 1.0 + t + 3.0 * t * t + 15.0 * t ** 3 + 105.0 * t ** 4
    return inv * s
