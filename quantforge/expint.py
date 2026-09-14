"""Exponential integrals E1, Ei, and the generalized E_n.

The exponential integral appears in radiative transfer, hydrology, and asymptotics. ``E1(x)
= integral_x^inf e^-t / t dt`` for ``x > 0``; ``Ei(x) = -integral_{-x}^inf e^-t / t dt``
(principal value) is its companion on the real line, with ``Ei(x) = -E1(-x)`` off the branch.
``E_n(x) = integral_1^inf e^{-x t} / t^n dt`` generalizes to integer order. Computed by the
convergent power series for small ``x`` and the continued fraction for large ``x``. Pure
standard library.
"""

import math

_EULER = 0.5772156649015328606


def e1(x):
    """Exponential integral ``E1(x) = integral_x^inf e^-t / t dt`` for ``x > 0``."""
    if x <= 0:
        raise ValueError("E1 requires x > 0")
    if x < 1.0:
        # power series: E1(x) = -gamma - ln x - sum_{k>=1} (-1)^k x^k / (k k!)
        s = 0.0
        term = 1.0
        for k in range(1, 100):
            term *= -x / k
            delta = term / k
            s += delta
            if abs(delta) < 1e-18 * abs(s):
                break
        return -_EULER - math.log(x) - s
    # continued fraction (Lentz) for x >= 1
    tiny = 1e-300
    b = x + 1.0
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 200):
        a = -i * i
        b += 2.0
        d = 1.0 / (b + a * d) if (b + a * d) != 0 else 1.0 / tiny
        c = b + a / c if (b + a / c) != 0 else tiny
        delta = c * d
        h *= delta
        if abs(delta - 1.0) < 1e-16:
            break
    return math.exp(-x) * h


def ei(x):
    """Exponential integral ``Ei(x)`` (principal value) for real ``x != 0``."""
    if x == 0:
        raise ValueError("Ei is singular at 0")
    if x < 0:
        return -e1(-x)
    # series for Ei(x), x > 0: Ei(x) = gamma + ln x + sum_{k>=1} x^k / (k k!)
    if x < 40.0:
        s = 0.0
        term = 1.0
        for k in range(1, 200):
            term *= x / k
            delta = term / k
            s += delta
            if abs(delta) < 1e-18 * abs(s):
                break
        return _EULER + math.log(x) + s
    # asymptotic series for large x
    s = 1.0
    term = 1.0
    for k in range(1, 40):
        term *= k / x
        s += term
        if term < 1e-16:
            break
    return math.exp(x) / x * s


def en(n, x):
    """Generalized exponential integral ``E_n(x) = integral_1^inf e^{-x t} / t^n dt``.

    ``n >= 0`` integer, ``x >= 0`` (``x > 0`` when ``n <= 1``). Uses ``E_0 = e^-x / x`` and
    the upward recurrence ``E_{n}(x) = (e^-x - x E_{n-1}(x)) / (n - 1)`` seeded from ``E_1``.
    """
    if n < 0:
        raise ValueError("n must be non-negative")
    if x < 0:
        raise ValueError("x must be non-negative")
    if n == 0:
        if x == 0:
            raise ValueError("E_0 singular at 0")
        return math.exp(-x) / x
    if n == 1:
        return e1(x)
    if x == 0:
        return 1.0 / (n - 1)   # E_n(0) = 1/(n-1) for n > 1
    val = e1(x)
    for k in range(2, n + 1):
        val = (math.exp(-x) - x * val) / (k - 1)
    return val
