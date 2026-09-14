"""Vector algebra: dot, cross, norms, angles, projection, reflection.

The everyday operations on real vectors -- dot and (3-D) cross products, Euclidean norm
and normalization, the angle between two vectors, the projection of one onto another, and
reflection about a plane normal. ``dot``/``norm``/``angle_between``/``project`` work in
any dimension; ``cross`` is 3-D. Pure standard library.
"""

import math


def dot(a, b):
    """Dot product of two equal-length vectors."""
    if len(a) != len(b):
        raise ValueError("vectors must have equal length")
    return sum(a[i] * b[i] for i in range(len(a)))


def cross(a, b):
    """Cross product of two 3-vectors: a vector perpendicular to both."""
    if len(a) != 3 or len(b) != 3:
        raise ValueError("cross product is defined for 3-vectors")
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def norm(a):
    """Euclidean (L2) norm of a vector."""
    return math.sqrt(sum(c * c for c in a))


def normalize(a):
    """Unit vector in the direction of ``a``; raises on the zero vector."""
    n = norm(a)
    if n == 0.0:
        raise ValueError("cannot normalize the zero vector")
    return tuple(c / n for c in a)


def angle_between(a, b):
    """Angle between two non-zero vectors in radians, in ``[0, pi]``.

    Uses the numerically stable ``atan2(|a x b|, a . b)`` form in 3-D and the clamped
    ``acos`` of the normalized dot product otherwise.
    """
    na = norm(a)
    nb = norm(b)
    if na == 0.0 or nb == 0.0:
        raise ValueError("angle undefined for the zero vector")
    if len(a) == 3 and len(b) == 3:
        return math.atan2(norm(cross(a, b)), dot(a, b))
    c = dot(a, b) / (na * nb)
    return math.acos(max(-1.0, min(1.0, c)))


def vector_project(a, b):
    """Vector projection of ``a`` onto ``b`` (the component of ``a`` along ``b``)."""
    nb2 = dot(b, b)
    if nb2 == 0.0:
        raise ValueError("cannot project onto the zero vector")
    s = dot(a, b) / nb2
    return tuple(s * bi for bi in b)


def vector_reject(a, b):
    """Vector rejection of ``a`` from ``b``: the component of ``a`` perpendicular to ``b``."""
    p = vector_project(a, b)
    return tuple(a[i] - p[i] for i in range(len(a)))


def reflect(a, normal):
    """Reflect vector ``a`` about the plane with unit-normalizable ``normal``.

    ``a - 2 (a . n_hat) n_hat`` where ``n_hat`` is the normalized normal -- the standard
    mirror reflection (e.g. a ray bouncing off a surface).
    """
    n_hat = normalize(normal)
    d = 2.0 * dot(a, n_hat)
    return tuple(a[i] - d * n_hat[i] for i in range(len(a)))
