"""Geometry on the unit sphere: angular distance, slerp, mean direction.

Directions in 3-D live on the unit sphere, and the natural operations there are angular
rather than Euclidean: the angle between two directions, the great-circle (spherical-linear)
interpolation between them, and the mean direction of a set. These sit under attitude
estimation, spherical statistics, and camera/orientation work. Vectors are 3-tuples and are
normalized internally where a unit vector is required. Pure standard library.
"""

import math


def _norm(v):
    return math.sqrt(sum(c * c for c in v))


def _normalize(v):
    n = _norm(v)
    if n == 0:
        raise ValueError("cannot normalize the zero vector")
    return [c / n for c in v]


def angular_distance(u, v):
    """Angle in radians between two 3-vectors (their directions), in ``[0, pi]``.

    Uses ``atan2(|u x v|, u . v)``, which stays accurate for both nearly-parallel and
    nearly-opposite directions (unlike ``acos`` of the dot product).
    """
    a = _normalize(u)
    b = _normalize(v)
    dot = sum(a[i] * b[i] for i in range(3))
    cx = [a[1] * b[2] - a[2] * b[1],
          a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]]
    return math.atan2(_norm(cx), dot)


def slerp_vectors(u, v, t):
    """Spherical-linear interpolation between unit directions ``u`` and ``v`` at ``t`` in [0,1].

    Returns a unit vector on the great-circle arc from ``u`` (``t=0``) to ``v`` (``t=1``),
    at constant angular speed. Falls back to normalized linear interpolation when the two
    directions are nearly identical.
    """
    a = _normalize(u)
    b = _normalize(v)
    dot = max(-1.0, min(1.0, sum(a[i] * b[i] for i in range(3))))
    theta = math.acos(dot)
    if theta < 1e-9:
        # nearly parallel: lerp then renormalize
        lin = [a[i] + t * (b[i] - a[i]) for i in range(3)]
        return _normalize(lin)
    s = math.sin(theta)
    w1 = math.sin((1 - t) * theta) / s
    w2 = math.sin(t * theta) / s
    return [w1 * a[i] + w2 * b[i] for i in range(3)]


def spherical_centroid(vectors):
    """Mean direction of a set of 3-vectors: the normalized vector sum.

    Returns the unit vector minimizing the sum of squared chord distances (the resultant
    direction of directional statistics). Raises if the vectors sum to zero (no mean
    direction).
    """
    if not vectors:
        raise ValueError("need at least one vector")
    total = [0.0, 0.0, 0.0]
    for v in vectors:
        u = _normalize(v)
        for i in range(3):
            total[i] += u[i]
    if _norm(total) < 1e-15:
        raise ValueError("vectors cancel; mean direction is undefined")
    return _normalize(total)


def spherical_resultant_length(vectors):
    """Mean resultant length ``R`` in ``[0, 1]``: concentration of a set of directions.

    ``R = |sum unit(v)| / n``. ``R = 1`` means all directions coincide; ``R = 0`` means they
    are perfectly spread. The dispersion measure of spherical statistics.
    """
    if not vectors:
        raise ValueError("need at least one vector")
    total = [0.0, 0.0, 0.0]
    for v in vectors:
        u = _normalize(v)
        for i in range(3):
            total[i] += u[i]
    return _norm(total) / len(vectors)
