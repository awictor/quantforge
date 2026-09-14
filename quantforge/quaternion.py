"""Unit quaternions for 3-D rotation.

Quaternions represent rotations without the gimbal lock of Euler angles and interpolate
smoothly (slerp) where rotation matrices cannot. A quaternion is ``(w, x, y, z)``; a unit
quaternion encodes a rotation of ``2*acos(w)`` about the axis ``(x, y, z)``. Provides the
Hamilton product, normalization, conversion to/from axis-angle, vector rotation, and
spherical linear interpolation. Pure standard library.
"""

import math


def quat_multiply(a, b):
    """Hamilton product ``a * b`` of two quaternions ``(w, x, y, z)``.

    Composition of rotations: ``quat_multiply(a, b)`` applies ``b`` then ``a``. Not
    commutative.
    """
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return (
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    )


def quat_normalize(q):
    """Return ``q`` scaled to unit norm; raises on the zero quaternion."""
    n = math.sqrt(sum(c * c for c in q))
    if n == 0.0:
        raise ValueError("cannot normalize the zero quaternion")
    return tuple(c / n for c in q)


def quat_conjugate(q):
    """Conjugate ``(w, -x, -y, -z)`` -- the inverse rotation for a unit quaternion."""
    w, x, y, z = q
    return (w, -x, -y, -z)


def axis_angle_to_quat(axis, angle):
    """Unit quaternion for a rotation of ``angle`` radians about ``axis`` (a 3-vector).

    The axis is normalized internally. ``angle = 0`` gives the identity ``(1, 0, 0, 0)``.
    """
    ax, ay, az = axis
    n = math.sqrt(ax * ax + ay * ay + az * az)
    if n == 0.0:
        raise ValueError("axis must be non-zero")
    s = math.sin(angle / 2.0) / n
    return (math.cos(angle / 2.0), ax * s, ay * s, az * s)


def quat_to_axis_angle(q):
    """Recover ``(axis, angle)`` from a unit quaternion; ``angle`` in ``[0, pi]``.

    Returns a unit ``axis`` and the rotation angle in radians. For the identity the axis
    is arbitrary; ``(1, 0, 0)`` is returned.
    """
    w, x, y, z = quat_normalize(q)
    angle = 2.0 * math.acos(max(-1.0, min(1.0, w)))
    s = math.sqrt(1.0 - w * w)
    if s < 1e-12:
        return (1.0, 0.0, 0.0), 0.0
    return (x / s, y / s, z / s), angle


def rotate_vector(q, v):
    """Rotate a 3-vector ``v`` by unit quaternion ``q`` (``q v q*``)."""
    q = quat_normalize(q)
    vq = (0.0, v[0], v[1], v[2])
    w, x, y, z = quat_multiply(quat_multiply(q, vq), quat_conjugate(q))
    return (x, y, z)


def slerp(a, b, t):
    """Spherical linear interpolation between unit quaternions ``a`` and ``b`` at ``t``.

    ``t = 0`` returns ``a``, ``t = 1`` returns ``b``, and intermediate ``t`` traces the
    shortest constant-speed arc on the unit sphere -- the standard smooth rotation blend.
    Chooses the shorter path (negates ``b`` if the dot product is negative).
    """
    a = quat_normalize(a)
    b = quat_normalize(b)
    dot = sum(a[i] * b[i] for i in range(4))
    if dot < 0.0:
        b = tuple(-c for c in b)
        dot = -dot
    if dot > 0.9995:
        # Nearly parallel: linear interpolate and renormalize.
        result = tuple(a[i] + t * (b[i] - a[i]) for i in range(4))
        return quat_normalize(result)
    theta0 = math.acos(dot)
    theta = theta0 * t
    sin0 = math.sin(theta0)
    s_a = math.sin(theta0 - theta) / sin0
    s_b = math.sin(theta) / sin0
    return tuple(s_a * a[i] + s_b * b[i] for i in range(4))
