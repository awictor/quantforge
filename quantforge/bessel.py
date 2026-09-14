"""Bessel functions of the first and second kind: J0, J1, Jn, Y0, Y1.

The Bessel functions solve the radial part of the wave and diffusion equations in cylindrical
geometry -- vibrating drumheads, waveguides, heat in a rod. ``J_n`` (first kind) is finite at
the origin; ``Y_n`` (second kind) diverges there. J0/J1/Y0/Y1 use the Abramowitz & Stegun
polynomial approximations (small argument) and the amplitude-phase asymptotic form (large
argument); higher ``J_n`` come from stable downward recurrence (Miller's algorithm). Pure
standard library.
"""

import math


def bessel_j0(x):
    """Bessel function of the first kind, order 0."""
    ax = abs(x)
    if ax < 8.0:
        y = x * x
        num = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
              + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))))
        den = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
              + y * (59272.64853 + y * (267.8532712 + y * 1.0))))
        return num / den
    z = 8.0 / ax
    y = z * z
    xx = ax - 0.785398164
    p0 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
         + y * (-0.2073370639e-5 + y * 0.2093887211e-6)))
    p1 = -0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5
         + y * (0.7621095161e-6 + y * (-0.934935152e-7))))
    return math.sqrt(0.636619772 / ax) * (math.cos(xx) * p0 - z * math.sin(xx) * p1)


def bessel_j1(x):
    """Bessel function of the first kind, order 1."""
    ax = abs(x)
    if ax < 8.0:
        y = x * x
        num = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
              + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))))
        den = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
              + y * (99447.43394 + y * (376.9991397 + y * 1.0))))
        return num / den
    z = 8.0 / ax
    y = z * z
    xx = ax - 2.356194491
    p0 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
         + y * (0.2457520174e-5 + y * (-0.240337019e-6))))
    p1 = 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5
         + y * (-0.88228987e-6 + y * 0.105787412e-6)))
    ans = math.sqrt(0.636619772 / ax) * (math.cos(xx) * p0 - z * math.sin(xx) * p1)
    return ans if x >= 0 else -ans


def bessel_y0(x):
    """Bessel function of the second kind, order 0 (``x > 0``)."""
    if x <= 0:
        raise ValueError("Y0 requires x > 0")
    if x < 8.0:
        y = x * x
        num = -2957821389.0 + y * (7062834065.0 + y * (-512359803.6
              + y * (10879881.29 + y * (-86327.92757 + y * 228.4622733))))
        den = 40076544269.0 + y * (745249964.8 + y * (7189466.438
              + y * (47447.26470 + y * (226.1030244 + y * 1.0))))
        return num / den + 0.636619772 * bessel_j0(x) * math.log(x)
    z = 8.0 / x
    y = z * z
    xx = x - 0.785398164
    p0 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
         + y * (-0.2073370639e-5 + y * 0.2093887211e-6)))
    p1 = -0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5
         + y * (0.7621095161e-6 + y * (-0.934935152e-7))))
    return math.sqrt(0.636619772 / x) * (math.sin(xx) * p0 + z * math.cos(xx) * p1)


def bessel_y1(x):
    """Bessel function of the second kind, order 1 (``x > 0``)."""
    if x <= 0:
        raise ValueError("Y1 requires x > 0")
    if x < 8.0:
        y = x * x
        num = x * (-4.900604943e13 + y * (1.275274390e13 + y * (-5.153438139e11
              + y * (7.349264551e9 + y * (-4.237922726e7 + y * 8.511937935e4)))))
        den = 2.499580570e14 + y * (4.244419664e12 + y * (3.733650367e10
              + y * (2.245904002e8 + y * (1.020426050e6 + y * (3.549632885e3 + y)))))
        return num / den + 0.636619772 * (bessel_j1(x) * math.log(x) - 1.0 / x)
    z = 8.0 / x
    y = z * z
    xx = x - 2.356194491
    p0 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
         + y * (0.2457520174e-5 + y * (-0.240337019e-6))))
    p1 = 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5
         + y * (-0.88228987e-6 + y * 0.105787412e-6)))
    return math.sqrt(0.636619772 / x) * (math.sin(xx) * p0 + z * math.cos(xx) * p1)


def bessel_jn(n, x):
    """Bessel function of the first kind, integer order ``n`` (``n >= 0``).

    ``n = 0, 1`` dispatch to the direct approximations; higher orders use Miller's stable
    downward recurrence, normalized by ``J0`` (or the ``sum`` identity).
    """
    if n < 0:
        raise ValueError("n must be non-negative")
    if n == 0:
        return bessel_j0(x)
    if n == 1:
        return bessel_j1(x)
    if x == 0.0:
        return 0.0
    ax = abs(x)
    tox = 2.0 / ax
    # downward recurrence from a high even starting order (Miller's algorithm, NR bessj)
    m = 2 * ((n + int(math.sqrt(60.0 * n)) + int(ax)) // 2 + 10)
    bjp = 0.0
    bj = 1.0
    ans = 0.0
    total = 0.0
    jsum = False           # accumulate every other term for the normalization sum
    for j in range(m, 0, -1):
        bjm = j * tox * bj - bjp
        bjp = bj
        bj = bjm
        if abs(bj) > 1e10:
            bj *= 1e-10
            bjp *= 1e-10
            ans *= 1e-10
            total *= 1e-10
        if jsum:
            total += bj
        jsum = not jsum
        if j == n:
            ans = bjp
    total = 2.0 * total - bj      # normalization: J0 + 2(J2 + J4 + ...) = 1
    ans /= total
    if x < 0 and n % 2 == 1:
        ans = -ans
    return ans
