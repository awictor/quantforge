"""Cubic interpolation: natural cubic spline and monotone Hermite.

The natural cubic spline (second derivative zero at both ends) gives a smooth
C2 curve through the knots, while the Fritsch-Carlson monotone cubic Hermite
preserves the monotonicity of the input data (no overshoot) -- important for
discount-factor and survival curves that must stay monotone. Both interpolate
strictly within the knot range. Pure standard library.
"""

import bisect


def _validate_knots(xs, ys):
    if len(xs) != len(ys):
        raise ValueError("xs and ys must have equal length")
    if len(xs) < 2:
        raise ValueError("need at least two knots")
    for i in range(len(xs) - 1):
        if xs[i + 1] <= xs[i]:
            raise ValueError("xs must be strictly increasing")


def natural_cubic_spline(xs, ys):
    """Build a natural cubic spline interpolant, returning a callable ``f(x)``.

    Solves the tridiagonal system for the second derivatives with zero-curvature
    (natural) end conditions. The returned function evaluates the piecewise cubic
    and is exact at the knots, C2 in between. Clamps to the end segments outside
    ``[xs[0], xs[-1]]``.
    """
    _validate_knots(xs, ys)
    n = len(xs)
    h = [xs[i + 1] - xs[i] for i in range(n - 1)]
    # Tridiagonal system for second derivatives m (m[0]=m[n-1]=0).
    alpha = [0.0] * n
    for i in range(1, n - 1):
        alpha[i] = 3.0 * ((ys[i + 1] - ys[i]) / h[i]
                          - (ys[i] - ys[i - 1]) / h[i - 1])
    l = [1.0] + [0.0] * (n - 1)
    mu = [0.0] * n
    z = [0.0] * n
    for i in range(1, n - 1):
        l[i] = 2.0 * (xs[i + 1] - xs[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]
    l[n - 1] = 1.0
    m = [0.0] * n
    for j in range(n - 2, -1, -1):
        m[j] = z[j] - mu[j] * m[j + 1]
    xs_list = list(xs)

    def f(x):
        i = bisect.bisect_right(xs_list, x) - 1
        i = min(max(i, 0), n - 2)
        dx = x - xs[i]
        a = ys[i]
        b = (ys[i + 1] - ys[i]) / h[i] - h[i] * (2.0 * m[i] + m[i + 1]) / 3.0
        c = m[i]
        d = (m[i + 1] - m[i]) / (3.0 * h[i])
        return a + dx * (b + dx * (c + dx * d))

    return f


def monotone_cubic(xs, ys):
    """Fritsch-Carlson monotone cubic Hermite interpolant, returning ``f(x)``.

    Chooses the Hermite tangents so the interpolant preserves the monotonicity of
    the data: no overshoot between knots. Exact at the knots. The standard choice
    for discount-factor / survival curves that must not wiggle below/above the data.
    """
    _validate_knots(xs, ys)
    n = len(xs)
    h = [xs[i + 1] - xs[i] for i in range(n - 1)]
    delta = [(ys[i + 1] - ys[i]) / h[i] for i in range(n - 1)]
    # Initial tangents: average of neighbouring secants.
    tang = [0.0] * n
    tang[0] = delta[0]
    tang[n - 1] = delta[n - 2]
    for i in range(1, n - 1):
        if delta[i - 1] * delta[i] <= 0.0:
            tang[i] = 0.0   # local extremum -> flat tangent (no overshoot)
        else:
            tang[i] = (delta[i - 1] + delta[i]) / 2.0
    # Fritsch-Carlson limiter to enforce monotonicity.
    for i in range(n - 1):
        if delta[i] == 0.0:
            tang[i] = 0.0
            tang[i + 1] = 0.0
        else:
            a = tang[i] / delta[i]
            b = tang[i + 1] / delta[i]
            s = a * a + b * b
            if s > 9.0:
                t = 3.0 / (s ** 0.5)
                tang[i] = t * a * delta[i]
                tang[i + 1] = t * b * delta[i]
    xs_list = list(xs)

    def f(x):
        i = bisect.bisect_right(xs_list, x) - 1
        i = min(max(i, 0), n - 2)
        t = (x - xs[i]) / h[i]
        t2 = t * t
        t3 = t2 * t
        h00 = 2.0 * t3 - 3.0 * t2 + 1.0
        h10 = t3 - 2.0 * t2 + t
        h01 = -2.0 * t3 + 3.0 * t2
        h11 = t3 - t2
        return (h00 * ys[i] + h10 * h[i] * tang[i]
                + h01 * ys[i + 1] + h11 * h[i] * tang[i + 1])

    return f
