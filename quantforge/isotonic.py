"""Isotonic (monotone) regression by the pool-adjacent-violators algorithm.

Fit the monotone step function closest to the data in weighted least squares: find
``y_hat`` minimizing ``sum w_i (y_i - y_hat_i)^2`` subject to ``y_hat`` being
non-decreasing (or non-increasing). The pool-adjacent-violators algorithm (PAVA)
solves this exactly in a single ``O(n)`` sweep -- whenever two adjacent blocks
violate the order, it pools them into their weighted mean and back-propagates. Used
for probability calibration (isotonic calibration), dose-response curves, and any
fit where the response is known to move one way. Pure standard library.
"""


def isotonic_regression(y, weights=None, increasing=True):
    """Weighted isotonic regression by pool-adjacent-violators.

    Returns the fitted values ``y_hat`` (same length as ``y``) forming the monotone
    sequence closest to ``y`` in weighted least squares. ``weights`` defaults to all
    ones; ``increasing=False`` fits a non-increasing sequence. The fit is a step
    function: tied blocks share their common weighted mean.
    """
    n = len(y)
    if n == 0:
        return []
    if weights is None:
        weights = [1.0] * n
    if len(weights) != n:
        raise ValueError("weights must match y in length")
    if any(w <= 0 for w in weights):
        raise ValueError("weights must be positive")

    # Fit non-decreasing; for non-increasing, negate, fit, negate back.
    sign = 1.0 if increasing else -1.0
    vals = [sign * v for v in y]

    # Each active block: [weighted mean, total weight, count]. Merge left while the
    # new block's mean is below its left neighbour's (a monotonicity violation).
    means = []
    ws = []
    counts = []
    for i in range(n):
        m = vals[i]
        w = weights[i]
        c = 1
        while means and means[-1] >= m:
            pw = ws[-1]
            pm = means.pop()
            ws.pop()
            pc = counts.pop()
            m = (pm * pw + m * w) / (pw + w)
            w = pw + w
            c = pc + c
        means.append(m)
        ws.append(w)
        counts.append(c)

    # Expand blocks back to per-point fitted values.
    out = []
    for m, c in zip(means, counts):
        out.extend([sign * m] * c)
    return out


def isotonic_fit(x, y, weights=None, increasing=True):
    """Fit an isotonic step function of ``x`` and return a predictor callable.

    Sorts by ``x``, runs :func:`isotonic_regression` on the reordered response, and
    returns ``(y_hat, predict)`` where ``y_hat`` is the fit aligned to the *original*
    input order and ``predict(x_new)`` interpolates the monotone step fit at a new
    point (linear between fitted knots, clamped to the endpoints). Ties in ``x`` share
    a fitted value.
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have equal length")
    if n == 0:
        raise ValueError("need at least one point")
    if weights is None:
        weights = [1.0] * n
    order = sorted(range(n), key=lambda i: x[i])
    xs = [x[i] for i in order]
    ys = [y[i] for i in order]
    wsorted = [weights[i] for i in order]
    fitted_sorted = isotonic_regression(ys, wsorted, increasing)

    # Align back to original order.
    y_hat = [0.0] * n
    for pos, i in enumerate(order):
        y_hat[i] = fitted_sorted[pos]

    def predict(x_new):
        if x_new <= xs[0]:
            return fitted_sorted[0]
        if x_new >= xs[-1]:
            return fitted_sorted[-1]
        lo, hi = 0, n - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if xs[mid] <= x_new:
                lo = mid
            else:
                hi = mid
        x0, x1 = xs[lo], xs[hi]
        f0, f1 = fitted_sorted[lo], fitted_sorted[hi]
        if x1 == x0:
            return f0
        t = (x_new - x0) / (x1 - x0)
        return f0 + t * (f1 - f0)

    return y_hat, predict
