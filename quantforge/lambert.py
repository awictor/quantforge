"""Lambert W function: the inverse of ``w -> w e^w``.

Lambert W solves ``w e^w = x`` for ``w`` -- the function that untangles equations where the
unknown appears both inside and outside an exponential (compound growth, delay equations,
enzyme kinetics). It is multivalued for ``-1/e <= x < 0``: the principal branch ``W0``
(``w >= -1``) and the secondary branch ``W_{-1}`` (``w <= -1``). Both are computed by a good
initial guess refined with Halley's cubic iteration. Pure standard library.
"""

import math

_INV_E = 1.0 / math.e


def lambert_w0(x):
    """Principal branch ``W0(x)`` solving ``w e^w = x`` with ``w >= -1``. Defined for ``x >= -1/e``."""
    if x < -_INV_E - 1e-12:
        raise ValueError("W0 is real only for x >= -1/e")
    if x == 0.0:
        return 0.0
    # initial guess
    if x < 0:
        # near the branch point use the series in sqrt(2(ex+1))
        p = math.sqrt(2 * (math.e * x + 1))
        w = -1 + p - p * p / 3
    elif x < 3:
        w = 0.5 * x if x < 1 else math.log(x + 1)
    else:
        lx = math.log(x)
        w = lx - math.log(lx)
    return _halley(x, w)


def lambert_wm1(x):
    """Secondary branch ``W_{-1}(x)`` with ``w <= -1``. Defined for ``-1/e <= x < 0``."""
    if not (-_INV_E - 1e-12 <= x < 0):
        raise ValueError("W_{-1} is real only for -1/e <= x < 0")
    if x >= -1e-300:
        return float("-inf")
    # initial guess: asymptotic w ~ ln(-x) - ln(-ln(-x))
    l1 = math.log(-x)
    l2 = math.log(-l1) if l1 < 0 else 0.0
    if x < -0.1:
        # near the branch point, use the sqrt series (lower branch)
        p = -math.sqrt(2 * (math.e * x + 1))
        w = -1 + p - p * p / 3
    else:
        w = l1 - l2
    return _halley(x, w)


def _halley(x, w, max_iter=60):
    for _ in range(max_iter):
        ew = math.exp(w)
        f = w * ew - x
        # Halley step for f(w) = w e^w - x
        denom = ew * (w + 1) - (w + 2) * f / (2 * w + 2) if w != -1 else ew * (w + 1)
        if denom == 0:
            break
        w_next = w - f / denom
        if abs(w_next - w) <= 1e-16 * (1 + abs(w_next)):
            w = w_next
            break
        w = w_next
    return w
