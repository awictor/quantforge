"""Hurst exponent via rescaled-range (R/S) analysis.

The Hurst exponent ``H`` measures the long-memory / persistence of a series:

  * ``H = 0.5``  -- no memory (a random walk's increments; white noise).
  * ``H > 0.5``  -- persistent / trending: moves tend to continue.
  * ``H < 0.5``  -- anti-persistent / mean-reverting: moves tend to reverse.

Rescaled-range analysis estimates it from how the rescaled range ``R/S`` of the
series grows with the window length ``n``. For a self-similar process
``E[R/S] ~ c * n^H``, so ``H`` is the slope of ``log(R/S)`` against ``log(n)``:

  1. split the series into non-overlapping windows of length ``n``;
  2. in each window, cumulatively sum the demeaned data, take the range (max minus
     min) of that cumulative profile, and divide by the window's standard
     deviation -- that is the window's rescaled range;
  3. average over windows to get ``(R/S)(n)``;
  4. regress ``log(R/S)`` on ``log(n)`` over a range of window sizes; the slope is
     ``H``.

Pure standard library.
"""

import math


def _mean(x):
    return sum(x) / len(x)


def rescaled_range(window):
    """Rescaled range ``R/S`` of a single window.

    ``R`` is the range of the cumulative demeaned series; ``S`` is the window's
    population standard deviation. Returns 0 when the window is constant.
    """
    n = len(window)
    if n < 2:
        raise ValueError("window needs at least 2 points")
    m = _mean(window)
    dev = [v - m for v in window]
    # Cumulative deviation profile.
    cum = []
    running = 0.0
    for d in dev:
        running += d
        cum.append(running)
    r = max(cum) - min(cum)
    s = math.sqrt(sum(d * d for d in dev) / n)
    if s == 0.0:
        return 0.0
    return r / s


def _ols_slope(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    den = sum((xs[i] - mx) ** 2 for i in range(n))
    if den == 0.0:
        raise ValueError("degenerate window sizes")
    return num / den


def hurst_exponent(x, min_window=8, max_window=None):
    """Estimate the Hurst exponent of ``x`` by rescaled-range analysis.

    Parameters
    ----------
    x : sequence of float
        The series (levels for a walk, or increments -- interpret ``H`` relative
        to what you feed in).
    min_window : int
        Smallest window length used in the log-log regression (>= 2).
    max_window : int, optional
        Largest window length. Defaults to ``len(x) // 2``.

    Returns
    -------
    float
        The Hurst exponent, the slope of ``log(R/S)`` versus ``log(n)`` over
        dyadic window sizes. ~0.5 for a memoryless series, >0.5 persistent,
        <0.5 mean-reverting.
    """
    n_total = len(x)
    if n_total < 2 * min_window:
        raise ValueError("series too short for the requested min_window")
    if min_window < 2:
        raise ValueError("min_window must be >= 2")
    if max_window is None:
        max_window = n_total // 2
    if max_window <= min_window:
        raise ValueError("max_window must exceed min_window")

    logs_n = []
    logs_rs = []
    n = min_window
    while n <= max_window:
        # Non-overlapping windows of length n.
        n_windows = n_total // n
        rs_vals = []
        for w in range(n_windows):
            window = x[w * n:(w + 1) * n]
            rs = rescaled_range(window)
            if rs > 0.0:
                rs_vals.append(rs)
        if rs_vals:
            logs_n.append(math.log(n))
            logs_rs.append(math.log(sum(rs_vals) / len(rs_vals)))
        n *= 2

    if len(logs_n) < 2:
        raise ValueError("not enough window sizes to estimate a slope")
    return _ols_slope(logs_n, logs_rs)
