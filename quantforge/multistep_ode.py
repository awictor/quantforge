"""Adams-Bashforth-Moulton predictor-corrector ODE integrator (4th order).

Where a Runge-Kutta step evaluates the derivative several times per step, a linear multistep
method reuses derivatives from previous steps: the 4-step Adams-Bashforth formula *predicts*
the next value, and the Adams-Moulton formula *corrects* it -- one new derivative evaluation
per step after the corrector (PECE), versus four for RK4. The first three steps are
bootstrapped with RK4. Solves ``y' = f(t, y)`` for scalar or vector ``y``. Pure standard
library.
"""


def _add(y, dy, h):
    if isinstance(y, (int, float)):
        return y + h * dy
    return [yi + h * di for yi, di in zip(y, dy)]


def _rk4_step(f, t, y, h):
    k1 = f(t, y)
    k2 = f(t + h / 2, _add(y, k1, h / 2))
    k3 = f(t + h / 2, _add(y, k2, h / 2))
    k4 = f(t + h, _add(y, k3, h))
    if isinstance(y, (int, float)):
        return y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return [yi + h / 6 * (a + 2 * b + 2 * c + d)
            for yi, a, b, c, d in zip(y, k1, k2, k3, k4)]


def _combine(y, terms):
    # y + sum of (coef, deriv) contributions
    if isinstance(y, (int, float)):
        return y + sum(c * d for c, d in terms)
    return [yi + sum(c * d[k] for c, d in terms) for k, yi in enumerate(y)]


def adams_bashforth_moulton(f, y0, t0, t1, n_steps):
    """Integrate ``y' = f(t, y)`` from ``t0`` to ``t1`` in ``n_steps`` using ABM4 (PECE).

    ``y0`` is scalar or a list (vector system). Returns ``(ts, ys)``: the ``n_steps + 1``
    time points and the solution at each. The first three steps use RK4 to build history.
    """
    if n_steps < 1:
        raise ValueError("n_steps must be at least 1")
    h = (t1 - t0) / n_steps
    ts = [t0 + i * h for i in range(n_steps + 1)]
    ys = [y0]
    fs = [f(t0, y0)]
    # bootstrap with RK4 until we have 4 derivative history points
    boot = min(3, n_steps)
    for i in range(boot):
        y_next = _rk4_step(f, ts[i], ys[i], h)
        ys.append(y_next)
        fs.append(f(ts[i + 1], y_next))
    # Adams-Bashforth-Moulton for the rest
    for i in range(boot, n_steps):
        # predictor: AB4 using f[i], f[i-1], f[i-2], f[i-3]
        yp = _combine(ys[i], [
            (h * 55 / 24, fs[i]),
            (-h * 59 / 24, fs[i - 1]),
            (h * 37 / 24, fs[i - 2]),
            (-h * 9 / 24, fs[i - 3]),
        ])
        fp = f(ts[i + 1], yp)
        # corrector: AM4 using the predicted derivative and recent history
        yc = _combine(ys[i], [
            (h * 9 / 24, fp),
            (h * 19 / 24, fs[i]),
            (-h * 5 / 24, fs[i - 1]),
            (h * 1 / 24, fs[i - 2]),
        ])
        ys.append(yc)
        fs.append(f(ts[i + 1], yc))
    return ts, ys
