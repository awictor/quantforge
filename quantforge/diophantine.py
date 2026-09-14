"""Diophantine equations: linear solutions, sqrt continued fractions, Pell.

Integer-solution tools. ``linear_diophantine`` solves ``a x + b y = c`` (all solutions
are one particular solution plus multiples of the homogeneous step).
``sqrt_continued_fraction`` gives the eventually-periodic continued-fraction expansion of
``sqrt(n)`` for non-square ``n``. ``pell_fundamental`` finds the smallest positive
solution of ``x^2 - n y^2 = 1`` from that expansion. Exact integer arithmetic. Pure
standard library.
"""

import math

from .modular import extended_gcd


def linear_diophantine(a, b, c):
    """Solve ``a x + b y = c`` in integers, or return ``None`` if unsolvable.

    Returns ``(x0, y0, dx, dy)``: one particular solution ``(x0, y0)`` and the step
    ``(dx, dy)`` so that ``(x0 + k*dx, y0 + k*dy)`` is a solution for every integer ``k``.
    Solvable iff ``gcd(a, b)`` divides ``c``. ``a`` and ``b`` must not both be zero.
    """
    if a == 0 and b == 0:
        raise ValueError("a and b cannot both be zero")
    g, x, y = extended_gcd(a, b)
    if c % g != 0:
        return None
    factor = c // g
    x0 = x * factor
    y0 = y * factor
    dx = b // g
    dy = -a // g
    return (x0, y0, dx, dy)


def sqrt_continued_fraction(n):
    """Continued-fraction expansion of ``sqrt(n)`` for a non-square ``n``: ``(a0, period)``.

    ``sqrt(n) = [a0; a1, a2, ...]`` is eventually periodic; returns the integer part
    ``a0`` and the repeating block ``period`` (a list). Raises if ``n`` is a perfect
    square or negative.
    """
    if n < 0:
        raise ValueError("n must be non-negative")
    a0 = math.isqrt(n)
    if a0 * a0 == n:
        raise ValueError("n must not be a perfect square")
    period = []
    m, d, a = 0, 1, a0
    while True:
        m = d * a - m
        d = (n - m * m) // d
        a = (a0 + m) // d
        period.append(a)
        if a == 2 * a0:              # standard termination for sqrt CFs
            break
    return a0, period


def pell_fundamental(n):
    """Smallest positive ``(x, y)`` solving Pell's equation ``x^2 - n y^2 = 1``.

    Built from the convergents of ``sqrt(n)``'s continued fraction. ``n`` must be a
    non-square positive integer.
    """
    if n <= 0:
        raise ValueError("n must be positive")
    a0 = math.isqrt(n)
    if a0 * a0 == n:
        raise ValueError("n must not be a perfect square")
    _, period = sqrt_continued_fraction(n)
    # Build convergents until x^2 - n y^2 == 1.
    # CF terms: a0, then the period repeated.
    def terms():
        yield a0
        while True:
            for t in period:
                yield t
    h_prev, h = 1, a0
    k_prev, k = 0, 1
    gen = terms()
    next(gen)                        # consume a0 (already in h)
    if h * h - n * k * k == 1:
        return (h, k)
    while True:
        a = next(gen)
        h_prev, h = h, a * h + h_prev
        k_prev, k = k, a * k + k_prev
        if h * h - n * k * k == 1:
            return (h, k)
