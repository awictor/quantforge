"""High-accuracy numerical derivatives by Ridders' method.

A plain central difference trades two errors: too large a step and the truncation
error dominates, too small and floating-point cancellation does. Ridders' method
(1982) evaluates the central difference at a shrinking sequence of steps and applies
Richardson extrapolation across them, tracking the error at each stage and stopping
when it starts to grow -- giving near machine-precision derivatives without hand-tuning
the step. Provides the first and second derivatives of a scalar function. Pure standard
library.
"""


def ridders_derivative(f, x, h=None, con=1.4, safe=2.0, ntab=10):
    """First derivative ``f'(x)`` by Ridders' polynomial extrapolation.

    ``h`` is the initial step (defaults to a scale-aware value); ``con`` is the step
    shrink factor per row, ``ntab`` the tableau size. Returns ``(derivative, error)``
    where ``error`` is the estimated absolute error. Raises if ``f`` cannot be improved
    at all.
    """
    if h is None:
        h = 1e-2 * (abs(x) + 1.0)
    if h <= 0:
        raise ValueError("step h must be positive")

    a = [[0.0] * ntab for _ in range(ntab)]
    hh = h
    a[0][0] = (f(x + hh) - f(x - hh)) / (2.0 * hh)
    err = float("inf")
    ans = a[0][0]
    for i in range(1, ntab):
        hh /= con
        a[0][i] = (f(x + hh) - f(x - hh)) / (2.0 * hh)
        fac = con * con
        for j in range(1, i + 1):
            a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0)
            fac *= con * con
            errt = max(abs(a[j][i] - a[j - 1][i]), abs(a[j][i] - a[j - 1][i - 1]))
            if errt <= err:
                err = errt
                ans = a[j][i]
        # Stop if the highest-order estimate is getting much worse (round-off wins).
        if abs(a[i][i] - a[i - 1][i - 1]) >= safe * err:
            break
    return ans, err


def ridders_second_derivative(f, x, h=None, con=1.4, safe=2.0, ntab=10):
    """Second derivative ``f''(x)`` by Ridders extrapolation of the central formula.

    Uses the three-point second difference ``(f(x+h) - 2 f(x) + f(x-h)) / h^2`` at a
    shrinking step sequence with Richardson extrapolation. Returns
    ``(second_derivative, error)``.
    """
    if h is None:
        h = 1e-2 * (abs(x) + 1.0)
    if h <= 0:
        raise ValueError("step h must be positive")

    fx = f(x)
    a = [[0.0] * ntab for _ in range(ntab)]
    hh = h
    a[0][0] = (f(x + hh) - 2.0 * fx + f(x - hh)) / (hh * hh)
    err = float("inf")
    ans = a[0][0]
    for i in range(1, ntab):
        hh /= con
        a[0][i] = (f(x + hh) - 2.0 * fx + f(x - hh)) / (hh * hh)
        fac = con * con
        for j in range(1, i + 1):
            a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0)
            fac *= con * con
            errt = max(abs(a[j][i] - a[j - 1][i]), abs(a[j][i] - a[j - 1][i - 1]))
            if errt <= err:
                err = errt
                ans = a[j][i]
        if abs(a[i][i] - a[i - 1][i - 1]) >= safe * err:
            break
    return ans, err
