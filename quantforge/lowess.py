"""LOWESS: locally-weighted scatterplot smoothing.

A nonparametric regression that fits a low-order polynomial in a moving neighbourhood,
weighting nearby points more via a tricube kernel. Unlike a global fit it follows
arbitrary curvature; unlike a spline it needs no knot placement, just a smoothing span.
Optional robustifying iterations (Cleveland 1979) down-weight points with large
residuals, so a few outliers don't distort the local fits. Returns the smoothed value
at each input abscissa. Pure standard library.
"""


def _tricube(u):
    a = abs(u)
    return (1.0 - a ** 3) ** 3 if a < 1.0 else 0.0


def _wls_local(xs, ys, weights, x0):
    """Weighted linear fit through a neighbourhood, evaluated at ``x0``."""
    sw = sum(weights)
    if sw == 0:
        return sum(ys) / len(ys)
    swx = sum(weights[i] * xs[i] for i in range(len(xs)))
    swy = sum(weights[i] * ys[i] for i in range(len(xs)))
    swxx = sum(weights[i] * xs[i] * xs[i] for i in range(len(xs)))
    swxy = sum(weights[i] * xs[i] * ys[i] for i in range(len(xs)))
    denom = sw * swxx - swx * swx
    if abs(denom) < 1e-300:
        return swy / sw
    slope = (sw * swxy - swx * swy) / denom
    intercept = (swy - slope * swx) / sw
    return slope * x0 + intercept


def _median(v):
    s = sorted(v)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else 0.5 * (s[mid - 1] + s[mid])


def lowess(x, y, frac=0.3, iterations=3):
    """LOWESS smooth of ``(x, y)``; returns the fitted value at each ``x``.

    ``frac`` is the span -- the fraction of points in each local neighbourhood (larger
    = smoother). ``iterations`` robustifying passes (Cleveland) down-weight outliers by
    a bisquare of their residuals; ``iterations=1`` disables robustifying. Points need
    not be sorted. Returns a list aligned to the input order.
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have equal length")
    if n < 2:
        raise ValueError("need at least 2 points")
    if not (0.0 < frac <= 1.0):
        raise ValueError("frac must be in (0, 1]")
    r = max(2, int(frac * n))              # neighbourhood size

    order = sorted(range(n), key=lambda i: x[i])
    xs = [x[i] for i in order]
    ys = [y[i] for i in order]

    robust = [1.0] * n
    fitted = [0.0] * n
    for _ in range(max(1, iterations)):
        for k in range(n):
            x0 = xs[k]
            # Distances to all points; pick the r nearest.
            dists = sorted(range(n), key=lambda j: abs(xs[j] - x0))[:r]
            dmax = max(abs(xs[j] - x0) for j in dists) or 1.0
            wx = [xs[j] for j in dists]
            wy = [ys[j] for j in dists]
            w = [_tricube((xs[j] - x0) / dmax) * robust[j] for j in dists]
            fitted[k] = _wls_local(wx, wy, w, x0)
        # Update robustness weights from residuals (bisquare).
        resid = [ys[k] - fitted[k] for k in range(n)]
        s = _median([abs(rr) for rr in resid])
        if s <= 0:
            # More than half the points fit exactly; anything with a nonzero
            # residual is an outlier and gets zero weight (median-scale is degenerate).
            if all(rr == 0.0 for rr in resid):
                break
            robust = [1.0 if rr == 0.0 else 0.0 for rr in resid]
            continue
        robust = []
        for rr in resid:
            u = rr / (6.0 * s)
            robust.append((1.0 - u * u) ** 2 if abs(u) < 1.0 else 0.0)

    # Map fitted values back to the original input order.
    out = [0.0] * n
    for pos, i in enumerate(order):
        out[i] = fitted[pos]
    return out
