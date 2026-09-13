"""Thiele's continued-fraction rational interpolation.

Where Neville/Newton fit a *polynomial* through points, Thiele's method fits a
*rational* function -- a ratio of polynomials -- expressed as a continued fraction

    R(x) = y0 + (x-x0) / ( rho1 + (x-x1) / ( rho2 + (x-x2) / ( ... ) ) ),

whose coefficients rho are the *inverse (reciprocal) differences* of the data. A
rational interpolant models poles and flattening tails that a polynomial of the same
node count cannot, and it still passes through every node exactly. This builds the
reciprocal-difference coefficients once and evaluates the continued fraction at any
point. Pure standard library.
"""


def thiele_coefficients(xs, ys):
    """Reciprocal-difference coefficients for Thiele's continued fraction.

    Returns the list ``rho`` whose entries are the leading inverse differences used by
    :func:`thiele_eval`. Length equals ``len(xs)``. Raises on duplicate nodes or a
    degenerate (zero) difference that stalls the recursion.
    """
    n = len(xs)
    if n != len(ys):
        raise ValueError("xs and ys must have equal length")
    if n < 2:
        raise ValueError("need at least 2 points")

    # rho[i][j] triangular table of inverse differences; keep full table then read
    # the diagonal rho[0][j].
    rho = [[0.0] * n for _ in range(n)]
    for i in range(n):
        rho[i][0] = ys[i]
    for i in range(n - 1):
        denom = ys[i] - ys[i + 1]
        if denom == 0.0:
            raise ValueError("degenerate difference (equal y at distinct x)")
        rho[i][1] = (xs[i] - xs[i + 1]) / denom
    for j in range(2, n):
        for i in range(n - j):
            denom = rho[i][j - 1] - rho[i + 1][j - 1]
            if denom == 0.0:
                raise ValueError("degenerate reciprocal difference; try reordering nodes")
            rho[i][j] = (xs[i] - xs[i + j]) / denom + rho[i + 1][j - 2]
    return [rho[0][j] for j in range(n)]


def thiele_eval(xs, coeffs, x):
    """Evaluate Thiele's continued fraction with reciprocal-difference ``coeffs``.

    ``coeffs`` come from :func:`thiele_coefficients`. Evaluated from the bottom up; the
    ``rho[0]`` and ``rho[1]`` entries seed the fraction. Returns the rational
    interpolant's value at ``x``.
    """
    n = len(coeffs)
    if n < 2:
        raise ValueError("need at least 2 coefficients")
    # Continued fraction: R = c0 + (x-x0)/(c1 + (x-x1)/((c2-c0) + ... )), where the
    # denominator at level k is (c_k - c_{k-2}) for k >= 2 and c_1 for k = 1. Evaluate
    # bottom-up.
    a = coeffs[n - 1] - (coeffs[n - 3] if (n - 1) >= 2 else 0.0)
    for k in range(n - 2, 0, -1):
        base = coeffs[k] - (coeffs[k - 2] if k >= 2 else 0.0)
        a = base + (x - xs[k]) / a
    return coeffs[0] + (x - xs[0]) / a


def thiele_interpolate(xs, ys, x):
    """Convenience: build Thiele coefficients and evaluate at ``x`` in one call."""
    return thiele_eval(xs, thiele_coefficients(xs, ys), x)
