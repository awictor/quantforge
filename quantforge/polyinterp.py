"""Polynomial interpolation: Neville's algorithm and Newton divided differences.

Given ``n`` points, there is a unique degree-``(n-1)`` polynomial through them. Two
classic ways to work with it:

  * ``neville`` -- evaluates that polynomial at a single point ``x`` by Neville's
    recursive tableau, returning both the value and an error estimate (the last
    correction applied), useful for Richardson-style extrapolation to a limit.
  * ``divided_differences`` / ``newton_polynomial`` -- build the Newton form once
    (``O(n^2)``) and then evaluate it cheaply at many points; adding a new point only
    appends one coefficient.

Pure standard library.
"""


def neville(xs, ys, x):
    """Evaluate the interpolating polynomial at ``x`` by Neville's algorithm.

    Returns ``(value, error_estimate)`` where the error estimate is the magnitude of
    the last correction -- a practical indicator of interpolation accuracy and the
    basis of Richardson extrapolation (interpolating to ``x = 0`` in the step size).
    """
    n = len(xs)
    if n != len(ys):
        raise ValueError("xs and ys must have equal length")
    if n == 0:
        raise ValueError("need at least one point")
    c = list(ys)
    d = list(ys)
    # Start from the nearest node for a well-conditioned path.
    ns = min(range(n), key=lambda i: abs(x - xs[i]))
    value = ys[ns]
    ns -= 1
    err = 0.0
    for m in range(1, n):
        for i in range(n - m):
            ho = xs[i] - x
            hp = xs[i + m] - x
            w = c[i + 1] - d[i]
            den = ho - hp
            if den == 0.0:
                raise ValueError("duplicate x values")
            den = w / den
            d[i] = hp * den
            c[i] = ho * den
        if 2 * (ns + 1) < (n - m):
            err = c[ns + 1]
        else:
            err = d[ns]
            ns -= 1
        value += err
    return value, abs(err)


def divided_differences(xs, ys):
    """Newton divided-difference coefficients for the points ``(xs, ys)``.

    Returns the list ``[f[x0], f[x0,x1], ...]`` -- the leading coefficients of the
    Newton form, computed in ``O(n^2)``. Feed these to :func:`newton_polynomial`.
    """
    n = len(xs)
    if n != len(ys):
        raise ValueError("xs and ys must have equal length")
    if n == 0:
        raise ValueError("need at least one point")
    coef = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            den = xs[i] - xs[i - j]
            if den == 0.0:
                raise ValueError("duplicate x values")
            coef[i] = (coef[i] - coef[i - 1]) / den
    return coef


def newton_polynomial(xs, coef, x):
    """Evaluate the Newton form with divided-difference ``coef`` at ``x`` (Horner)."""
    n = len(coef)
    p = coef[-1]
    for k in range(n - 2, -1, -1):
        p = p * (x - xs[k]) + coef[k]
    return p
