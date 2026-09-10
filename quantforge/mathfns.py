"""Special functions used by the pricing engine, implemented from scratch.

We avoid SciPy/NumPy so the library has zero third-party dependencies. The
normal CDF is built on a high-accuracy erf approximation; the tests pin its
error against reference values.
"""

import math

SQRT_2 = math.sqrt(2.0)
SQRT_2PI = math.sqrt(2.0 * math.pi)
INV_SQRT_2PI = 1.0 / SQRT_2PI


def norm_pdf(x: float) -> float:
    """Standard normal probability density function."""
    return INV_SQRT_2PI * math.exp(-0.5 * x * x)


def norm_cdf(x: float) -> float:
    """Standard normal cumulative distribution function.

    Uses the C standard-library ``math.erf`` (available since Python 3.2),
    which is correct to full double precision. ``Phi(x) = 0.5 * erfc(-x/sqrt2)``.
    """
    return 0.5 * math.erfc(-x / SQRT_2)


def norm_ppf(p: float) -> float:
    """Inverse standard normal CDF (quantile function).

    Acklam's rational approximation refined by one Halley step against the
    exact ``norm_cdf``. Relative error after refinement is < 1e-15 across the
    open interval (0, 1).
    """
    if not (0.0 < p < 1.0):
        if p == 0.0:
            return -math.inf
        if p == 1.0:
            return math.inf
        raise ValueError("norm_ppf requires 0 <= p <= 1")

    # Coefficients for Acklam's algorithm.
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]

    plow = 0.02425
    phigh = 1.0 - plow

    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    elif p <= phigh:
        q = p - 0.5
        r = q * q
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
    else:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)

    # One Halley refinement step for full double precision.
    e = norm_cdf(x) - p
    u = e * SQRT_2PI * math.exp(0.5 * x * x)
    x = x - u / (1.0 + 0.5 * x * u)
    return x
