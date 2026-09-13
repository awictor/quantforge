"""Bivariate and trivariate standard-normal CDFs.

The cumulative probabilities of correlated standard normals, needed to price
multi-asset and compound options in closed form:

  * ``bivariate_normal_cdf`` -- ``P(X1 <= a, X2 <= b)`` with correlation ``rho``, by the
    Drezner-Wesolowsky single-integral form (the routine behind the library's American
    and compound-option models, exposed here).
  * ``trivariate_normal_cdf`` -- ``P(X1 <= a, X2 <= b, X3 <= c)`` for a 3x3 correlation
    matrix, by reducing to a one-dimensional integral of the bivariate CDF (Genz).

Pure standard library.
"""

import math

from .american import _bivariate_normal
from .mathfns import norm_cdf
from .quadrature import _GL


def bivariate_normal_cdf(a, b, rho):
    """Standard bivariate normal CDF ``P(X1 <= a, X2 <= b; corr=rho)``.

    ``rho`` in ``[-1, 1]``. Accurate to ~1e-7 across the usable correlation range.
    """
    if not (-1.0 <= rho <= 1.0):
        raise ValueError("rho must be in [-1, 1]")
    return _bivariate_normal(a, b, rho)


def trivariate_normal_cdf(a, b, c, r12, r13, r23, n=24):
    """Standard trivariate normal CDF ``P(X1<=a, X2<=b, X3<=c)``.

    ``r12, r13, r23`` are the pairwise correlations of a valid 3x3 correlation matrix.
    Uses the reduction ``P3 = P(X3<=c) * ...`` via a 1-D integral over the third
    variable of a conditional bivariate CDF (Genz). ``n`` sets the quadrature panels.
    Falls back to the product/independent forms when correlations vanish.
    """
    for r in (r12, r13, r23):
        if not (-1.0 <= r <= 1.0):
            raise ValueError("correlations must be in [-1, 1]")
    if abs(r12) < 1e-12 and abs(r13) < 1e-12 and abs(r23) < 1e-12:
        return norm_cdf(a) * norm_cdf(b) * norm_cdf(c)

    # Integrate over x3 = z in (-inf, c]: contribution phi(z) * P(X1<=a, X2<=b | X3=z).
    # Conditional on X3=z: mean1 = r13 z, var1 = 1-r13^2; mean2 = r23 z, var2 = 1-r23^2;
    # conditional corr rho_c = (r12 - r13 r23) / sqrt((1-r13^2)(1-r23^2)).
    s1 = math.sqrt(max(1.0 - r13 * r13, 1e-300))
    s2 = math.sqrt(max(1.0 - r23 * r23, 1e-300))
    denom = s1 * s2
    rho_c = (r12 - r13 * r23) / denom if denom > 0 else 0.0
    rho_c = max(-1.0, min(1.0, rho_c))

    # Composite Simpson over a truncated lower tail up to c.
    lo = -8.0
    hi = c
    if hi <= lo:
        return 0.0
    m = max(2, n)
    if m % 2:
        m += 1
    h = (hi - lo) / m
    inv_sqrt_2pi = 1.0 / math.sqrt(2.0 * math.pi)

    def integrand(z):
        phi = inv_sqrt_2pi * math.exp(-0.5 * z * z)
        a_c = (a - r13 * z) / s1
        b_c = (b - r23 * z) / s2
        return phi * _bivariate_normal(a_c, b_c, rho_c)

    total = integrand(lo) + integrand(hi)
    for i in range(1, m):
        total += (4.0 if i % 2 else 2.0) * integrand(lo + i * h)
    val = total * h / 3.0
    return max(0.0, min(1.0, val))
