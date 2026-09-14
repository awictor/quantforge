"""CORDIC: sine, cosine, atan2, and magnitude by shift-and-add only.

CORDIC (COordinate Rotation DIgital Computer) evaluates trigonometric and vector functions
using only additions, subtractions, and bit shifts -- no multiplies -- by rotating a vector
through a fixed sequence of ever-smaller angles ``atan(2^-i)``. *Rotation* mode drives the
angle to a target to get ``(cos, sin)``; *vectoring* mode drives the y-component to zero to
get ``atan2`` and the magnitude. It is how fixed-point hardware and early calculators
computed these functions. Here it is a readable reference, cross-checked against ``math``.
Pure standard library.
"""

import math

_ITERS = 40
# precomputed rotation angles atan(2^-i) and the cumulative gain K
_ANGLES = [math.atan(2.0 ** -i) for i in range(_ITERS)]
_K = 1.0
for _a in range(_ITERS):
    _K *= 1.0 / math.sqrt(1.0 + 4.0 ** -_a)


def cordic_sincos(theta):
    """Return ``(cos(theta), sin(theta))`` by CORDIC rotation mode.

    ``theta`` is reduced into ``[-pi/2, pi/2]`` (the algorithm's convergence range) using the
    identities for the outer quadrants. Accurate to ~1e-10 with 40 iterations.
    """
    # reduce to [-pi, pi]
    t = math.fmod(theta, 2 * math.pi)
    if t > math.pi:
        t -= 2 * math.pi
    elif t < -math.pi:
        t += 2 * math.pi
    # fold into [-pi/2, pi/2], tracking a sign flip for the outer quadrants
    flip = 1.0
    if t > math.pi / 2:
        t -= math.pi
        flip = -1.0
    elif t < -math.pi / 2:
        t += math.pi
        flip = -1.0
    x, y, z = _K, 0.0, t
    for i in range(_ITERS):
        d = 1.0 if z >= 0 else -1.0
        x, y = x - d * y * (2.0 ** -i), y + d * x * (2.0 ** -i)
        z -= d * _ANGLES[i]
    return (flip * x, flip * y)


def cordic_atan2(y, x):
    """Return ``atan2(y, x)`` in ``(-pi, pi]`` by CORDIC vectoring mode."""
    if x == 0 and y == 0:
        return 0.0
    # vectoring converges for x > 0; handle other quadrants by reflection
    add = 0.0
    if x < 0:
        if y >= 0:
            add = math.pi
        else:
            add = -math.pi
        x, y = -x, -y
    vx, vy, vz = x, y, 0.0
    for i in range(_ITERS):
        d = -1.0 if vy >= 0 else 1.0
        vx, vy = vx - d * vy * (2.0 ** -i), vy + d * vx * (2.0 ** -i)
        vz -= d * _ANGLES[i]
    return vz + add


def cordic_hypot(x, y):
    """Return ``sqrt(x^2 + y^2)`` by CORDIC vectoring mode (no square root)."""
    if x == 0 and y == 0:
        return 0.0
    vx, vy = abs(x), abs(y)
    for i in range(_ITERS):
        d = -1.0 if vy >= 0 else 1.0
        vx, vy = vx - d * vy * (2.0 ** -i), vy + d * vx * (2.0 ** -i)
    return vx * _K
