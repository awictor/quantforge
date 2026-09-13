"""Special functions: incomplete gamma and beta, digamma, inverse error function.

The regularized incomplete gamma and beta functions and their relatives underpin
the gamma, chi-square, Poisson, Student-t, F and binomial distributions used across
the library. This module gives clean, public, well-tested implementations so those
distribution routines need not each re-derive them. All pure standard library on
top of :func:`math.lgamma`.

Conventions:

- ``gammainc(a, x)`` is the regularized lower incomplete gamma
  ``P(a, x) = gamma(a, x) / Gamma(a)``, and ``gammaincc(a, x) = 1 - P = Q(a, x)``.
- ``betainc(a, b, x)`` is the regularized incomplete beta ``I_x(a, b)``.
- ``digamma(x)`` is ``psi(x) = d/dx ln Gamma(x)``.
- ``erfinv(y)`` inverts the error function on ``(-1, 1)``.
"""

import math

_MAXIT = 300
_EPS = 1e-15
_TINY = 1e-300


def gammainc(a, x):
    """Regularized lower incomplete gamma ``P(a, x) = gamma(a, x) / Gamma(a)``.

    Series for ``x < a + 1`` and the complement of the Lentz continued fraction
    otherwise (Numerical Recipes). ``P(a, 0) = 0`` and ``P(a, x) -> 1`` as
    ``x -> inf``. Requires ``a > 0`` and ``x >= 0``.
    """
    if a <= 0.0:
        raise ValueError("a must be positive")
    if x < 0.0:
        raise ValueError("x must be non-negative")
    if x == 0.0:
        return 0.0
    gln = math.lgamma(a)
    if x < a + 1.0:
        ap = a
        term = 1.0 / a
        total = term
        for _ in range(_MAXIT):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * _EPS:
                break
        return total * math.exp(-x + a * math.log(x) - gln)
    return 1.0 - _gammaincc_cf(a, x, gln)


def _gammaincc_cf(a, x, gln):
    """Upper regularized gamma ``Q(a, x)`` by the Lentz continued fraction."""
    b = x + 1.0 - a
    c = 1.0 / _TINY
    d = 1.0 / b
    h = d
    for i in range(1, _MAXIT):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < _TINY:
            d = _TINY
        c = b + an / c
        if abs(c) < _TINY:
            c = _TINY
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < _EPS:
            break
    return math.exp(-x + a * math.log(x) - gln) * h


def gammaincc(a, x):
    """Regularized upper incomplete gamma ``Q(a, x) = 1 - P(a, x)``.

    Uses the continued fraction directly for ``x >= a + 1`` (where it converges
    fast) and the series complement otherwise, so ``gammainc(a, x) + gammaincc(a, x)
    == 1`` to machine precision.
    """
    if a <= 0.0:
        raise ValueError("a must be positive")
    if x < 0.0:
        raise ValueError("x must be non-negative")
    if x == 0.0:
        return 1.0
    if x < a + 1.0:
        return 1.0 - gammainc(a, x)
    return _gammaincc_cf(a, x, math.lgamma(a))


def betainc(a, b, x):
    """Regularized incomplete beta ``I_x(a, b)`` via the Lentz continued fraction.

    ``I_0 = 0``, ``I_1 = 1``, and the symmetry ``I_x(a, b) = 1 - I_{1-x}(b, a)`` is
    used where the fraction converges slowly. Requires ``a, b > 0`` and
    ``0 <= x <= 1``.
    """
    if a <= 0.0 or b <= 0.0:
        raise ValueError("a and b must be positive")
    if not (0.0 <= x <= 1.0):
        raise ValueError("x must be in [0, 1]")
    if x == 0.0:
        return 0.0
    if x == 1.0:
        return 1.0
    # Evaluate the fraction where it converges fast; swap arguments otherwise.
    if x < (a + 1.0) / (a + b + 2.0):
        return _betacf_result(a, b, x)
    return 1.0 - _betacf_result(b, a, 1.0 - x)


def _betacf_result(a, b, x):
    """``I_x(a, b)`` from the Lentz continued fraction (fast-converging branch)."""
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(math.log(x) * a + math.log(1.0 - x) * b - lbeta) / a
    f, c, d = 1.0, 1.0, 0.0
    for i in range(0, 2 * _MAXIT):
        m = i // 2
        if i == 0:
            num = 1.0
        elif i % 2 == 0:
            num = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
        else:
            num = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1.0 + num * d
        if abs(d) < _TINY:
            d = _TINY
        d = 1.0 / d
        c = 1.0 + num / c
        if abs(c) < _TINY:
            c = _TINY
        f *= c * d
        if abs(1.0 - c * d) < _EPS:
            break
    return front * (f - 1.0)


def digamma(x):
    """Digamma ``psi(x) = d/dx ln Gamma(x)`` for ``x > 0``.

    Recurses up to ``x >= 6`` with ``psi(x) = psi(x+1) - 1/x``, then applies the
    asymptotic (Stirling) series. Accurate to ~1e-12 for positive arguments.
    """
    if x <= 0.0:
        raise ValueError("x must be positive")
    result = 0.0
    while x < 6.0:
        result -= 1.0 / x
        x += 1.0
    inv = 1.0 / x
    inv2 = inv * inv
    # psi(x) ~ ln x - 1/(2x) - sum Bernoulli terms.
    result += (math.log(x) - 0.5 * inv
               - inv2 * (1.0 / 12.0
                         - inv2 * (1.0 / 120.0
                                   - inv2 * (1.0 / 252.0 - inv2 / 240.0))))
    return result


def erfinv(y):
    """Inverse error function on ``(-1, 1)``.

    A rational approximation (Giles) seeds a Newton-Halley refinement against
    :func:`math.erf`, giving full double precision. ``erfinv(0) = 0`` and
    ``erfinv(erf(x)) = x``.
    """
    if not (-1.0 < y < 1.0):
        if y == 1.0:
            return math.inf
        if y == -1.0:
            return -math.inf
        raise ValueError("y must be in (-1, 1)")
    if y == 0.0:
        return 0.0
    w = -math.log((1.0 - y) * (1.0 + y))
    if w < 5.0:
        w -= 2.5
        p = 2.81022636e-08
        for c in (3.43273939e-07, -3.5233877e-06, -4.39150654e-06,
                  0.00021858087, -0.00125372503, -0.00417768164,
                  0.246640727, 1.50140941):
            p = p * w + c
    else:
        w = math.sqrt(w) - 3.0
        p = -0.000200214257
        for c in (0.000100950558, 0.00134934322, -0.00367342844,
                  0.00573950773, -0.0076224613, 0.00943887047,
                  1.00167406, 2.83297682):
            p = p * w + c
    x = p * y
    # Two Newton steps against erf to reach machine precision.
    two_over_sqrt_pi = 2.0 / math.sqrt(math.pi)
    for _ in range(2):
        err = math.erf(x) - y
        x -= err / (two_over_sqrt_pi * math.exp(-x * x))
    return x
