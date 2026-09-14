"""Conversions between the three common 3-D rotation representations.

A rotation can be written as a unit quaternion ``(w, x, y, z)``, a 3x3 orthogonal matrix, or
a triple of Euler angles. Each form suits a different task -- quaternions interpolate and
compose without gimbal lock, matrices apply directly to vectors, Euler angles read out as
yaw/pitch/roll -- so converting between them is routine. Euler angles here use the aerospace
Z-Y-X (yaw, pitch, roll) intrinsic convention. Quaternions match
:mod:`quantforge.quaternion` (``(w, x, y, z)``, scalar first). Pure standard library.
"""

import math


def quat_to_matrix(q):
    """Convert a unit quaternion ``(w, x, y, z)`` to a 3x3 rotation matrix (row lists)."""
    w, x, y, z = q
    n = math.sqrt(w * w + x * x + y * y + z * z)
    if n == 0:
        raise ValueError("zero quaternion has no rotation")
    w, x, y, z = w / n, x / n, y / n, z / n
    return [
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ]


def matrix_to_quat(m):
    """Convert a 3x3 rotation matrix to a unit quaternion ``(w, x, y, z)`` (Shepperd's method)."""
    trace = m[0][0] + m[1][1] + m[2][2]
    if trace > 0:
        s = 0.5 / math.sqrt(trace + 1.0)
        w = 0.25 / s
        x = (m[2][1] - m[1][2]) * s
        y = (m[0][2] - m[2][0]) * s
        z = (m[1][0] - m[0][1]) * s
    elif m[0][0] > m[1][1] and m[0][0] > m[2][2]:
        s = 2.0 * math.sqrt(1.0 + m[0][0] - m[1][1] - m[2][2])
        w = (m[2][1] - m[1][2]) / s
        x = 0.25 * s
        y = (m[0][1] + m[1][0]) / s
        z = (m[0][2] + m[2][0]) / s
    elif m[1][1] > m[2][2]:
        s = 2.0 * math.sqrt(1.0 + m[1][1] - m[0][0] - m[2][2])
        w = (m[0][2] - m[2][0]) / s
        x = (m[0][1] + m[1][0]) / s
        y = 0.25 * s
        z = (m[1][2] + m[2][1]) / s
    else:
        s = 2.0 * math.sqrt(1.0 + m[2][2] - m[0][0] - m[1][1])
        w = (m[1][0] - m[0][1]) / s
        x = (m[0][2] + m[2][0]) / s
        y = (m[1][2] + m[2][1]) / s
        z = 0.25 * s
    return _canonical((w, x, y, z))


def euler_to_quat(yaw, pitch, roll):
    """Z-Y-X intrinsic Euler angles (yaw, pitch, roll, radians) to a quaternion ``(w,x,y,z)``."""
    cy, sy = math.cos(yaw / 2), math.sin(yaw / 2)
    cp, sp = math.cos(pitch / 2), math.sin(pitch / 2)
    cr, sr = math.cos(roll / 2), math.sin(roll / 2)
    w = cr * cp * cy + sr * sp * sy
    x = sr * cp * cy - cr * sp * sy
    y = cr * sp * cy + sr * cp * sy
    z = cr * cp * sy - sr * sp * cy
    return (w, x, y, z)


def quat_to_euler(q):
    """Quaternion ``(w, x, y, z)`` to Z-Y-X Euler angles ``(yaw, pitch, roll)`` in radians."""
    w, x, y, z = q
    n = math.sqrt(w * w + x * x + y * y + z * z)
    if n == 0:
        raise ValueError("zero quaternion has no rotation")
    w, x, y, z = w / n, x / n, y / n, z / n
    # roll (x-axis)
    roll = math.atan2(2 * (w * x + y * z), 1 - 2 * (x * x + y * y))
    # pitch (y-axis), clamped for the gimbal-lock pole
    sinp = 2 * (w * y - z * x)
    sinp = max(-1.0, min(1.0, sinp))
    pitch = math.asin(sinp)
    # yaw (z-axis)
    yaw = math.atan2(2 * (w * z + x * y), 1 - 2 * (y * y + z * z))
    return (yaw, pitch, roll)


def euler_to_matrix(yaw, pitch, roll):
    """Z-Y-X Euler angles directly to a 3x3 rotation matrix."""
    return quat_to_matrix(euler_to_quat(yaw, pitch, roll))


def matrix_to_euler(m):
    """3x3 rotation matrix to Z-Y-X Euler angles ``(yaw, pitch, roll)``."""
    return quat_to_euler(matrix_to_quat(m))


def _canonical(q):
    # fix the sign so w >= 0 (q and -q are the same rotation) for stable comparison
    if q[0] < 0:
        return tuple(-c for c in q)
    return tuple(q)
