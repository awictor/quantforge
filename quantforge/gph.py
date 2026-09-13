"""Geweke-Porter-Hudak (GPH) estimator of the long-memory parameter d.

A fractionally-integrated series has a spectral density that behaves like
``f(lambda) ~ C * lambda^{-2d}`` near the origin, so ``log f`` is linear in
``log(4 sin^2(lambda/2))`` with slope ``-d``. The GPH estimator regresses the log
periodogram on that regressor over the lowest ``m`` Fourier frequencies:

    log I(lambda_j) = c - d * log(4 sin^2(lambda_j / 2)) + error,

and reads ``d`` off the slope. The OLS standard error uses the known error variance
``pi^2 / 6`` of the log-periodogram. Pairs with :mod:`quantforge.fracdiff`: estimate
``d`` here, then difference by it. Pure standard library.
"""

import math

from .spectral import periodogram


def gph_estimate(x, m=None):
    """Estimate the fractional-integration order ``d`` by the GPH regression.

    ``x`` is the series; ``m`` is the number of low frequencies used (defaults to
    ``floor(sqrt(n))``, the common bandwidth). Returns a dict with ``d`` (the memory
    parameter), ``std_error`` (asymptotic, from the ``pi^2/6`` log-periodogram
    variance), ``m`` and ``n``. ``d ~ 0`` indicates short memory, ``d > 0`` long
    memory, ``d < 0`` anti-persistence.
    """
    n = len(x)
    if n < 8:
        raise ValueError("need at least 8 observations")
    if m is None:
        m = int(math.sqrt(n))
    if m < 2 or m >= n // 2:
        raise ValueError("m must be in [2, n/2)")

    freqs, power = periodogram(x)
    # Use frequencies j = 1 .. m (skip the DC term at index 0).
    xs = []
    ys = []
    for j in range(1, m + 1):
        lam = 2.0 * math.pi * freqs[j]              # angular frequency
        reg = math.log(4.0 * math.sin(lam / 2.0) ** 2)
        if power[j] <= 0.0:
            continue
        xs.append(reg)
        ys.append(math.log(power[j]))
    if len(xs) < 2:
        raise ValueError("not enough positive periodogram ordinates")

    k = len(xs)
    mx = sum(xs) / k
    my = sum(ys) / k
    sxx = sum((xi - mx) ** 2 for xi in xs)
    sxy = sum((xs[i] - mx) * (ys[i] - my) for i in range(k))
    slope = sxy / sxx
    d = -slope
    # Asymptotic standard error: Var(log I) = pi^2/6, so Var(d_hat) = (pi^2/6)/Sxx.
    std_error = math.sqrt((math.pi ** 2 / 6.0) / sxx)
    return {"d": d, "std_error": std_error, "m": m, "n": n}


def fractional_integrate(noise, d):
    """Fractionally integrate white ``noise`` by order ``d``: apply ``(1 - L)^{-d}``.

    The inverse of :func:`quantforge.fracdiff.fractional_difference`; useful for
    generating ARFIMA(0, d, 0) test series. Uses the binomial weights of
    ``(1 - L)^{-d}``, ``w_0 = 1``, ``w_k = w_{k-1} (k - 1 + d) / k``.
    """
    n = len(noise)
    if n == 0:
        raise ValueError("noise must be non-empty")
    w = [1.0]
    for k in range(1, n):
        w.append(w[-1] * (k - 1 + d) / k)
    out = []
    for t in range(n):
        out.append(sum(w[k] * noise[t - k] for k in range(t + 1)))
    return out
