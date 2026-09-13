"""Ordinary differential equation solvers: RK4 and adaptive RK45.

Integrate an initial-value problem ``y'(t) = f(t, y)`` from ``t0`` to ``t1``:

  * ``rk4`` -- the classic fourth-order Runge-Kutta with a fixed step. Simple and
    accurate (global error ``O(h^4)``) when the solution is smooth and no step control
    is needed.
  * ``rk45`` -- the Dormand-Prince adaptive method: a 5th-order step with an embedded
    4th-order estimate gives a local error estimate that resizes the step to hold a
    tolerance, taking small steps only where the solution moves fast.

``y`` may be a scalar or a vector (list); ``f`` returns the same shape. Pure standard
library.
"""


def _add(y, k, s):
    return [y[i] + s * k[i] for i in range(len(y))]


def _as_vec(y):
    return list(y) if isinstance(y, (list, tuple)) else [y]


def rk4(f, t0, y0, t1, n=100):
    """Fixed-step RK4 from ``t0`` to ``t1`` in ``n`` steps.

    ``f(t, y)`` returns the derivative (scalar or list matching ``y0``). Returns
    ``(ts, ys)``: the ``n+1`` time points and the state at each (each state a list).
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    vec = isinstance(y0, (list, tuple))
    y = _as_vec(y0)
    d = len(y)
    h = (t1 - t0) / n

    def fv(t, yv):
        r = f(t, yv if vec else yv[0])
        return _as_vec(r)

    ts = [t0]
    ys = [list(y)]
    t = t0
    for _ in range(n):
        k1 = fv(t, y)
        k2 = fv(t + h / 2, _add(y, k1, h / 2))
        k3 = fv(t + h / 2, _add(y, k2, h / 2))
        k4 = fv(t + h, _add(y, k3, h))
        y = [y[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(d)]
        t += h
        ts.append(t)
        ys.append(list(y))
    return ts, ys


# Dormand-Prince coefficients.
_C = [0.0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0, 1.0]
_A = [
    [],
    [1 / 5],
    [3 / 40, 9 / 40],
    [44 / 45, -56 / 15, 32 / 9],
    [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
    [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
    [35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84],
]
_B5 = [35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0.0]
_B4 = [5179 / 57600, 0.0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40]


def rk45(f, t0, y0, t1, tol=1e-8, h0=None, max_steps=100000):
    """Adaptive Dormand-Prince (RK45) integration from ``t0`` to ``t1``.

    Controls the step to keep the estimated local error near ``tol``. Returns
    ``(ts, ys)`` at the accepted steps (non-uniform). ``h0`` is the initial step
    (defaults to a fraction of the interval).
    """
    vec = isinstance(y0, (list, tuple))
    y = _as_vec(y0)
    d = len(y)

    def fv(t, yv):
        return _as_vec(f(t, yv if vec else yv[0]))

    h = (t1 - t0) / 100.0 if h0 is None else h0
    t = t0
    ts = [t0]
    ys = [list(y)]
    steps = 0
    while t < t1 - 1e-15 and steps < max_steps:
        steps += 1
        if t + h > t1:
            h = t1 - t
        ks = []
        for i in range(7):
            yi = list(y)
            for j in range(i):
                yi = _add(yi, ks[j], h * _A[i][j])
            ks.append(fv(t + _C[i] * h, yi))
        y5 = [y[m] + h * sum(_B5[i] * ks[i][m] for i in range(7)) for m in range(d)]
        y4 = [y[m] + h * sum(_B4[i] * ks[i][m] for i in range(7)) for m in range(d)]
        err = max(abs(y5[m] - y4[m]) for m in range(d))
        if err <= tol or h <= 1e-14:
            t += h
            y = y5
            ts.append(t)
            ys.append(list(y))
        # Step-size update (safety factor 0.9, order-5 exponent).
        if err > 0:
            h *= min(4.0, max(0.1, 0.9 * (tol / err) ** 0.2))
        else:
            h *= 2.0
    return ts, ys
