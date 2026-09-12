"""One-dimensional minimizers: golden-section and Brent's method.

For a scalar unimodal function on a bracketing interval:

  * ``golden_section_min`` -- shrinks the bracket by the golden ratio each step;
    robust, derivative-free, linear convergence.
  * ``brent_min`` -- Brent's method, combining a parabolic-interpolation step with
    a golden-section safeguard for superlinear convergence when the function is
    smooth.

Both return ``(x_min, f_min)``. These are the line-search / 1-D calibration
building blocks (e.g. implied-parameter fitting). Pure standard library.
"""

import math

_INV_PHI = (math.sqrt(5.0) - 1.0) / 2.0          # 1/phi ~ 0.618
_INV_PHI2 = (3.0 - math.sqrt(5.0)) / 2.0         # 1/phi^2


def golden_section_min(f, lo, hi, tol=1e-10, max_iter=200):
    """Minimize a unimodal ``f`` on ``[lo, hi]`` by golden-section search."""
    if lo >= hi:
        raise ValueError("require lo < hi")
    a, b = lo, hi
    c = b - _INV_PHI * (b - a)
    d = a + _INV_PHI * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(max_iter):
        if abs(b - a) < tol:
            break
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - _INV_PHI * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + _INV_PHI * (b - a)
            fd = f(d)
    xm = 0.5 * (a + b)
    return xm, f(xm)


def brent_min(f, lo, hi, tol=1e-10, max_iter=200):
    """Minimize a scalar ``f`` on ``[lo, hi]`` by Brent's method.

    Parabolic interpolation with a golden-section fallback. Returns
    ``(x_min, f_min)``.
    """
    if lo >= hi:
        raise ValueError("require lo < hi")
    a, b = lo, hi
    x = w = v = a + _INV_PHI2 * (b - a)
    fx = fw = fv = f(x)
    e = 0.0
    d = 0.0
    for _ in range(max_iter):
        m = 0.5 * (a + b)
        tol1 = tol * abs(x) + 1e-12
        tol2 = 2.0 * tol1
        if abs(x - m) <= tol2 - 0.5 * (b - a):
            break
        use_golden = True
        if abs(e) > tol1:
            # Fit a parabola through (x, fx), (w, fw), (v, fv).
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2.0 * (q - r)
            if q > 0.0:
                p = -p
            q = abs(q)
            etemp = e
            e = d
            if not (abs(p) >= abs(0.5 * q * etemp) or p <= q * (a - x)
                    or p >= q * (b - x)):
                d = p / q
                u = x + d
                if u - a < tol2 or b - u < tol2:
                    d = tol1 if x < m else -tol1
                use_golden = False
        if use_golden:
            e = (b - x) if x < m else (a - x)
            d = _INV_PHI2 * e
        u = x + d if abs(d) >= tol1 else x + (tol1 if d > 0 else -tol1)
        fu = f(u)
        if fu <= fx:
            if u < x:
                b = x
            else:
                a = x
            v, fv, w, fw, x, fx = w, fw, x, fx, u, fu
        else:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v, fv, w, fw = w, fw, u, fu
            elif fu <= fv or v == x or v == w:
                v, fv = u, fu
    return x, fx
