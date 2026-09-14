"""Catmull-Rom splines: a C1 interpolating spline through a set of points.

Unlike a Bezier curve (which only touches its endpoints), a Catmull-Rom spline passes
*through* every control point, with each segment's tangents set from the neighbouring
points -- so it interpolates while staying smooth (C1). The ``alpha`` parameter selects the
parametrization: ``0`` uniform, ``0.5`` centripetal (no cusps or self-intersections, the
usual choice), ``1`` chordal. Points are scalars or equal-length coordinate tuples. Pure
standard library.
"""


def _sub(p, q):
    if isinstance(p, (int, float)):
        return p - q
    return tuple(a - b for a, b in zip(p, q))


def _dist(p, q):
    if isinstance(p, (int, float)):
        return abs(p - q)
    return sum((a - b) ** 2 for a, b in zip(p, q)) ** 0.5


def _segment(p0, p1, p2, p3, t, alpha):
    # non-uniform Catmull-Rom on the middle segment [p1, p2], parameter t in [0, 1]
    def tj(ti, pa, pb):
        d = _dist(pa, pb)
        return ti + d ** alpha if d > 0 else ti + 1e-12  # avoid coincident-point division

    t0 = 0.0
    t1 = tj(t0, p0, p1)
    t2 = tj(t1, p1, p2)
    t3 = tj(t2, p2, p3)
    tt = t1 + (t2 - t1) * t

    def lerp(a, b, ta, tb, x):
        w = (tb - x) / (tb - ta)
        if isinstance(a, (int, float)):
            return w * a + (1 - w) * b
        return tuple(w * ai + (1 - w) * bi for ai, bi in zip(a, b))

    a1 = lerp(p0, p1, t0, t1, tt)
    a2 = lerp(p1, p2, t1, t2, tt)
    a3 = lerp(p2, p3, t2, t3, tt)
    b1 = lerp(a1, a2, t0, t2, tt)
    b2 = lerp(a2, a3, t1, t3, tt)
    return lerp(b1, b2, t1, t2, tt)


def catmull_rom_point(points, seg, t, alpha=0.5):
    """Evaluate segment ``seg`` (between ``points[seg]`` and ``points[seg+1]``) at ``t`` in [0,1].

    Needs at least two points; endpoints are handled by duplicating the boundary point as the
    phantom neighbour. ``alpha``: 0 uniform, 0.5 centripetal, 1 chordal.
    """
    n = len(points)
    if n < 2:
        raise ValueError("need at least two points")
    if not 0 <= seg < n - 1:
        raise IndexError("segment index out of range")
    p1 = points[seg]
    p2 = points[seg + 1]
    p0 = points[seg - 1] if seg - 1 >= 0 else p1
    p3 = points[seg + 2] if seg + 2 < n else p2
    return _segment(p0, p1, p2, p3, t, alpha)


def catmull_rom_curve(points, samples_per_segment=16, alpha=0.5):
    """Sample the whole spline. Returns a flat list of points along all segments.

    Each segment contributes ``samples_per_segment`` points; the shared endpoints appear
    once (the start of each segment except the first is skipped). The curve passes through
    every input point.
    """
    n = len(points)
    if n < 2:
        raise ValueError("need at least two points")
    if samples_per_segment < 2:
        raise ValueError("samples_per_segment must be at least 2")
    out = [points[0]]
    for seg in range(n - 1):
        for k in range(1, samples_per_segment):
            t = k / (samples_per_segment - 1)
            out.append(catmull_rom_point(points, seg, t, alpha))
    return out
