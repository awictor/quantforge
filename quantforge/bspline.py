"""B-spline basis functions and curves (Cox-de Boor recursion, de Boor's algorithm).

A B-spline generalizes the Bezier curve: a piecewise-polynomial of degree ``p`` whose shape
is controlled by a knot vector and control points, with *local* control (moving one point
changes only a few spans). ``bspline_basis`` evaluates the Cox-de Boor basis functions,
``open_uniform_knots`` builds a clamped knot vector (so the curve touches its first and last
control points), and ``bspline_point`` evaluates the curve by de Boor's algorithm. Points
are scalars or equal-length coordinate tuples. Pure standard library.
"""


def open_uniform_knots(n_control, degree):
    """Clamped (open-uniform) knot vector for ``n_control`` points of the given ``degree``.

    Length ``n_control + degree + 1``: the first and last ``degree + 1`` knots are repeated
    (0 and 1), the interior knots equally spaced -- so the curve interpolates its endpoints.
    """
    if n_control < degree + 1:
        raise ValueError("need at least degree + 1 control points")
    p = degree
    m = n_control + p + 1
    knots = [0.0] * (p + 1)
    interior = n_control - p - 1
    for i in range(1, interior + 1):
        knots.append(i / (interior + 1))
    knots += [1.0] * (p + 1)
    assert len(knots) == m
    return knots


def bspline_basis(i, p, knots, t):
    """Cox-de Boor basis function ``N_{i,p}(t)`` for knot vector ``knots``.

    ``i`` is the basis index, ``p`` the degree. Uses the standard recursion with the
    ``0/0 = 0`` convention for repeated knots.
    """
    if p == 0:
        # half-open spans, but the very last span is closed so t == last knot evaluates
        if knots[i] <= t < knots[i + 1]:
            return 1.0
        if t == knots[-1] and knots[i] <= t <= knots[i + 1] and knots[i + 1] == knots[-1]:
            return 1.0
        return 0.0
    left = 0.0
    denom_l = knots[i + p] - knots[i]
    if denom_l > 0:
        left = (t - knots[i]) / denom_l * bspline_basis(i, p - 1, knots, t)
    right = 0.0
    denom_r = knots[i + p + 1] - knots[i + 1]
    if denom_r > 0:
        right = (knots[i + p + 1] - t) / denom_r * bspline_basis(i + 1, p - 1, knots, t)
    return left + right


def bspline_point(control, degree, t, knots=None):
    """Evaluate the B-spline curve at parameter ``t`` (basis-weighted control points).

    ``degree`` is the polynomial degree ``p``; ``knots`` defaults to a clamped open-uniform
    vector (so the curve interpolates the endpoints). ``t`` runs over the knot range
    ``[knots[p], knots[-p-1]]`` (``[0, 1]`` for the default knots).
    """
    n = len(control)
    if n < degree + 1:
        raise ValueError("need at least degree + 1 control points")
    if knots is None:
        knots = open_uniform_knots(n, degree)
    scalar = isinstance(control[0], (int, float))
    if scalar:
        acc = 0.0
    else:
        dim = len(control[0])
        acc = [0.0] * dim
    for i in range(n):
        b = bspline_basis(i, degree, knots, t)
        if b == 0.0:
            continue
        if scalar:
            acc += b * control[i]
        else:
            for k in range(dim):
                acc[k] += b * control[i][k]
    return acc if scalar else tuple(acc)


def bspline_curve(control, degree, samples, knots=None):
    """Sample the B-spline at ``samples`` equally spaced parameters over its knot range."""
    n = len(control)
    if samples < 2:
        raise ValueError("samples must be at least 2")
    if knots is None:
        knots = open_uniform_knots(n, degree)
    lo = knots[degree]
    hi = knots[len(knots) - degree - 1]
    out = []
    for j in range(samples):
        t = lo + (hi - lo) * j / (samples - 1)
        if j == samples - 1:
            t = hi
        out.append(bspline_point(control, degree, t, knots))
    return out
