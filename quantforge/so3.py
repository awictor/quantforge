"""SO(3) exponential and logarithm: rotation vectors and Rodrigues' formula.

A 3-D rotation can be written as a *rotation vector* ``omega`` -- direction is the axis,
magnitude is the angle. Rodrigues' formula exponentiates it to a rotation matrix,
``R = I + sin(t) K + (1 - cos t) K^2`` with ``K`` the skew-symmetric ``hat(omega/t)`` and
``t = |omega|``; the matrix logarithm inverts it. This is the exp/log map of the SO(3) Lie
group -- the natural way to interpolate, average, and take derivatives of rotations. Pure
standard library.
"""

import math


def hat(v):
    """Skew-symmetric matrix of a 3-vector: ``hat(v) w == cross(v, w)``."""
    x, y, z = v
    return [[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]]


def unhat(m):
    """Inverse of :func:`hat`: the 3-vector of a skew-symmetric matrix."""
    return (m[2][1], m[0][2], m[1][0])


def rodrigues(omega):
    """Exponential map: rotation vector ``omega`` -> 3x3 rotation matrix (Rodrigues' formula)."""
    x, y, z = omega
    theta = math.sqrt(x * x + y * y + z * z)
    if theta < 1e-12:
        # first-order expansion near zero: R ~ I + hat(omega)
        K = hat(omega)
        return [[(1.0 if i == j else 0.0) + K[i][j] for j in range(3)] for i in range(3)]
    kx, ky, kz = x / theta, y / theta, z / theta
    K = hat((kx, ky, kz))
    s = math.sin(theta)
    c = math.cos(theta)
    K2 = _matmul(K, K)
    R = [[(1.0 if i == j else 0.0) + s * K[i][j] + (1 - c) * K2[i][j]
          for j in range(3)] for i in range(3)]
    return R


def so3_log(R):
    """Logarithm map: rotation matrix ``R`` -> rotation vector ``omega`` (axis * angle)."""
    trace = R[0][0] + R[1][1] + R[2][2]
    cos_t = (trace - 1) / 2
    cos_t = max(-1.0, min(1.0, cos_t))
    theta = math.acos(cos_t)
    if theta < 1e-8:
        # near identity: omega ~ unhat of the skew part
        return (0.5 * (R[2][1] - R[1][2]),
                0.5 * (R[0][2] - R[2][0]),
                0.5 * (R[1][0] - R[0][1]))
    if abs(math.pi - theta) < 1e-6:
        # near pi: axis from the largest diagonal of (R + I)/2
        a = [(R[i][i] + 1) / 2 for i in range(3)]
        i = max(range(3), key=lambda k: a[i] if False else a[k])
        axis = [0.0, 0.0, 0.0]
        axis[i] = math.sqrt(max(0.0, a[i]))
        for j in range(3):
            if j != i:
                axis[j] = (R[i][j] + R[j][i]) / (4 * axis[i]) if axis[i] > 1e-12 else 0.0
        norm = math.sqrt(sum(c * c for c in axis))
        axis = [c / norm for c in axis] if norm > 0 else [1.0, 0.0, 0.0]
        return tuple(theta * c for c in axis)
    factor = theta / (2 * math.sin(theta))
    return (factor * (R[2][1] - R[1][2]),
            factor * (R[0][2] - R[2][0]),
            factor * (R[1][0] - R[0][1]))


def _matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
