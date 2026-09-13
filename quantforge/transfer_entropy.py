"""Information-theoretic dependence: mutual information and transfer entropy.

Estimated from histograms (fixed-bin discretization) of the series:

- ``mutual_information`` -- the reduction in uncertainty about one series from
  knowing another, ``I(X;Y) = sum p(x,y) log(p(x,y) / (p(x) p(y)))``; zero iff the
  two are independent, and symmetric in its arguments,
- ``transfer_entropy`` -- Schreiber's directed measure of how much ``X``'s past
  reduces uncertainty about ``Y``'s next value beyond ``Y``'s own past,
  ``TE_{X->Y} = sum p(y_{t+1}, y_t, x_t) log( p(y_{t+1}|y_t,x_t) / p(y_{t+1}|y_t) )``.
  Unlike mutual information it is directional, so ``TE_{X->Y}`` and ``TE_{Y->X}``
  differ when one series drives the other.

Logs are natural (nats). Pure standard library.
"""

import math


def _bin_index(value, lo, width, bins):
    if width <= 0.0:
        return 0
    k = int((value - lo) / width)
    return 0 if k < 0 else (bins - 1 if k >= bins else k)


def _discretize(series, bins):
    lo = min(series)
    hi = max(series)
    width = (hi - lo) / bins if hi > lo else 0.0
    return [_bin_index(v, lo, width, bins) for v in series]


def mutual_information(x, y, bins=8):
    """Mutual information ``I(X;Y)`` in nats, from a 2-D histogram.

    Non-negative, zero iff ``X`` and ``Y`` are independent, and symmetric:
    ``mutual_information(x, y) == mutual_information(y, x)``. ``bins`` sets the
    discretization resolution. Aligned series of at least two points.
    """
    n = len(x)
    if n < 2 or len(y) != n:
        raise ValueError("x and y must be equal-length with at least 2 points")
    if bins < 2:
        raise ValueError("bins must be at least 2")
    dx = _discretize(x, bins)
    dy = _discretize(y, bins)
    pxy = {}
    px = {}
    py = {}
    for i in range(n):
        pxy[(dx[i], dy[i])] = pxy.get((dx[i], dy[i]), 0) + 1
        px[dx[i]] = px.get(dx[i], 0) + 1
        py[dy[i]] = py.get(dy[i], 0) + 1
    mi = 0.0
    for (a, b), c in pxy.items():
        p_ab = c / n
        p_a = px[a] / n
        p_b = py[b] / n
        mi += p_ab * math.log(p_ab / (p_a * p_b))
    return max(0.0, mi)


def transfer_entropy(source, target, bins=8):
    """Transfer entropy ``TE_{source -> target}`` in nats (lag-1, Schreiber).

    Measures how much the source's present reduces uncertainty about the target's
    next value beyond the target's own present. Built from the joint histogram of
    ``(target_{t+1}, target_t, source_t)``. Directional: run it both ways to see
    which series leads. Non-negative; near zero when the source carries no extra
    information about the target's future. Aligned series of at least three points.
    """
    n = len(source)
    if n < 3 or len(target) != n:
        raise ValueError("source and target must be equal-length, >= 3 points")
    if bins < 2:
        raise ValueError("bins must be at least 2")
    ds = _discretize(source, bins)
    dt = _discretize(target, bins)

    # Counts over (y_{t+1}, y_t, x_t) for t = 0 .. n-2.
    p_yyx = {}   # (y1, y0, x0)
    p_yy = {}    # (y1, y0)
    p_yx = {}    # (y0, x0)
    p_y = {}     # (y0)
    m = n - 1
    for t in range(m):
        y1, y0, x0 = dt[t + 1], dt[t], ds[t]
        p_yyx[(y1, y0, x0)] = p_yyx.get((y1, y0, x0), 0) + 1
        p_yy[(y1, y0)] = p_yy.get((y1, y0), 0) + 1
        p_yx[(y0, x0)] = p_yx.get((y0, x0), 0) + 1
        p_y[y0] = p_y.get(y0, 0) + 1

    te = 0.0
    for (y1, y0, x0), c in p_yyx.items():
        p_joint = c / m
        # p(y1 | y0, x0) = p(y1,y0,x0) / p(y0,x0);  p(y1 | y0) = p(y1,y0) / p(y0)
        p_y1_given_yx = c / p_yx[(y0, x0)]
        p_y1_given_y = p_yy[(y1, y0)] / p_y[y0]
        te += p_joint * math.log(p_y1_given_yx / p_y1_given_y)
    return max(0.0, te)
