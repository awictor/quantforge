"""Detrended fluctuation analysis (DFA) of a time series.

DFA estimates a self-similarity scaling exponent that, unlike plain rescaled-range
analysis, is robust to slowly-varying trends and mild nonstationarity. The steps:

1. form the integrated profile ``Y_t = sum_{i<=t} (x_i - mean(x))``,
2. split ``Y`` into non-overlapping windows of length ``s``,
3. in each window fit and subtract a least-squares polynomial trend (order 1 here),
4. take the root-mean-square of the residuals as the fluctuation ``F(s)``,
5. repeat over a range of ``s`` and read the scaling exponent ``alpha`` off the slope
   of ``log F(s)`` against ``log s``.

For white noise ``alpha ~ 0.5``; ``alpha > 0.5`` is persistent (long-range
correlated), ``alpha < 0.5`` anti-persistent, and integrating a series adds 1 to its
exponent (a random walk gives ``alpha ~ 1.5``). Pure standard library.
"""

import math


def _detrend_rms(segment):
    """RMS of a segment after removing its least-squares linear trend."""
    n = len(segment)
    xs = list(range(n))
    mx = (n - 1) / 2.0
    my = sum(segment) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((xs[i] - mx) * (segment[i] - my) for i in range(n))
    slope = sxy / sxx if sxx > 0 else 0.0
    intercept = my - slope * mx
    ss = 0.0
    for i in range(n):
        resid = segment[i] - (slope * i + intercept)
        ss += resid * resid
    return math.sqrt(ss / n)


def dfa_fluctuations(x, scales=None):
    """Fluctuation function ``F(s)`` of DFA over a set of window sizes ``scales``.

    Returns ``(scales, fluctuations)``. ``scales`` defaults to a dyadic-ish grid
    between 4 and ``len(x)//4``. Windows that do not divide the series exactly drop
    the remainder. Requires at least 16 points.
    """
    n = len(x)
    if n < 16:
        raise ValueError("need at least 16 observations")
    mean = sum(x) / n
    profile = []
    acc = 0.0
    for xi in x:
        acc += xi - mean
        profile.append(acc)

    if scales is None:
        scales = []
        s = 4
        while s <= n // 4:
            scales.append(s)
            s = int(s * 1.5) if int(s * 1.5) > s else s + 1
    used = []
    flucts = []
    for s in scales:
        if s < 4 or s > n // 2:
            continue
        n_seg = n // s
        rms = []
        for w in range(n_seg):
            seg = profile[w * s:(w + 1) * s]
            rms.append(_detrend_rms(seg))
        f = math.sqrt(sum(r * r for r in rms) / len(rms))
        if f > 0.0:
            used.append(s)
            flucts.append(f)
    if len(used) < 2:
        raise ValueError("not enough valid scales")
    return used, flucts


def dfa_exponent(x, scales=None):
    """Detrended-fluctuation scaling exponent ``alpha``.

    The slope of ``log F(s)`` regressed on ``log s`` from :func:`dfa_fluctuations`.
    ``alpha ~ 0.5`` for white noise, ``> 0.5`` for persistent series, ``< 0.5`` for
    anti-persistent; it equals the Hurst exponent for a stationary long-memory
    series and exceeds it by 1 for the integrated (random-walk) version.
    """
    scales, flucts = dfa_fluctuations(x, scales)
    xs = [math.log(s) for s in scales]
    ys = [math.log(f) for f in flucts]
    k = len(xs)
    mx = sum(xs) / k
    my = sum(ys) / k
    sxx = sum((xi - mx) ** 2 for xi in xs)
    sxy = sum((xs[i] - mx) * (ys[i] - my) for i in range(k))
    return sxy / sxx
