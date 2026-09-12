"""General-purpose scalar root finders: bisection, Brent, Newton.

Robust one-dimensional root finding for the many implied-parameter solves in the
library (implied vol, yields, breakevens). Bisection is the guaranteed-convergent
fallback; Brent's method combines bisection, secant, and inverse quadratic
interpolation for fast robust convergence; Newton uses the derivative with a
bisection safeguard. Pure standard library.
"""


def bisection(f, lo, hi, tol=1e-12, max_iter=200):
    """Bisection root of ``f`` on ``[lo, hi]`` (requires a sign change).

    Halves the bracket until it is narrower than ``tol``. Guaranteed to converge
    for a continuous ``f`` with ``f(lo) f(hi) < 0``.
    """
    flo, fhi = f(lo), f(hi)
    if flo == 0.0:
        return lo
    if fhi == 0.0:
        return hi
    if flo * fhi > 0.0:
        raise ValueError("f(lo) and f(hi) must bracket a root (opposite signs)")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if abs(fmid) < tol or 0.5 * (hi - lo) < tol:
            return mid
        if flo * fmid < 0.0:
            hi = mid
        else:
            lo, flo = mid, fmid
    return 0.5 * (lo + hi)


def brent(f, lo, hi, tol=1e-12, max_iter=200):
    """Brent's method root of ``f`` on ``[lo, hi]`` (requires a sign change).

    Combines bisection with secant and inverse-quadratic interpolation for
    superlinear convergence while retaining bisection's guaranteed bracketing.
    """
    a, b = lo, hi
    fa, fb = f(a), f(b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0.0:
        raise ValueError("f(lo) and f(hi) must bracket a root (opposite signs)")
    if abs(fa) < abs(fb):
        a, b = b, a
        fa, fb = fb, fa
    c, fc = a, fa
    mflag = True
    d = a
    for _ in range(max_iter):
        if abs(b - a) < tol or fb == 0.0:
            return b
        if fa != fc and fb != fc:
            # Inverse quadratic interpolation.
            s = (a * fb * fc / ((fa - fb) * (fa - fc))
                 + b * fa * fc / ((fb - fa) * (fb - fc))
                 + c * fa * fb / ((fc - fa) * (fc - fb)))
        else:
            s = b - fb * (b - a) / (fb - fa)   # secant
        cond = (not (min((3 * a + b) / 4.0, b) <= s <= max((3 * a + b) / 4.0, b))
                or (mflag and abs(s - b) >= abs(b - c) / 2.0)
                or (not mflag and abs(s - b) >= abs(c - d) / 2.0)
                or (mflag and abs(b - c) < tol)
                or (not mflag and abs(c - d) < tol))
        if cond:
            s = 0.5 * (a + b)   # bisection step
            mflag = True
        else:
            mflag = False
        fs = f(s)
        d, c, fc = c, b, fb
        if fa * fs < 0.0:
            b, fb = s, fs
        else:
            a, fa = s, fs
        if abs(fa) < abs(fb):
            a, b = b, a
            fa, fb = fb, fa
    return b


def newton(f, fprime, x0, tol=1e-12, max_iter=100, lo=None, hi=None):
    """Newton's method with an optional bisection safeguard.

    Steps ``x -= f(x)/f'(x)``; if ``lo``/``hi`` bounds are given, a step leaving
    the bracket (or a zero derivative) falls back to a bisection step. Converges
    quadratically near a simple root when the derivative is well-behaved.
    """
    x = x0
    if lo is not None and hi is not None:
        flo = f(lo)
    for _ in range(max_iter):
        fx = f(x)
        if abs(fx) < tol:
            return x
        dfx = fprime(x)
        if dfx == 0.0:
            if lo is None or hi is None:
                raise ValueError("zero derivative and no bracket to fall back on")
            x = 0.5 * (lo + hi)
            continue
        x_new = x - fx / dfx
        if lo is not None and hi is not None:
            if not (lo < x_new < hi):
                x_new = 0.5 * (lo + hi)
            # Maintain the bracket.
            if flo * fx < 0.0:
                hi = x
            else:
                lo, flo = x, fx
        x = x_new
    return x
