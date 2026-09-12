"""Structural-break diagnostics: CUSUM of mean and the Chow test.

Financial relationships are not stationary forever -- a mean, a beta, or a
volatility regime can shift. Two classic detectors:

  * ``cusum_mean`` -- the standardized cumulative sum of deviations from the mean.
    Under a stable mean it wanders inside a confidence band that widens like the
    Brownian bridge; a persistent level shift makes the path drift out of the band.
  * ``chow_test`` -- an F test for whether the mean (or a regression) differs
    before and after a known break point. Large F rejects "no break".

Pure standard library.
"""

import math


def _mean(x):
    return sum(x) / len(x)


def _var(x, ddof=1):
    m = _mean(x)
    n = len(x)
    if n - ddof <= 0:
        raise ValueError("not enough points for the requested ddof")
    return sum((v - m) ** 2 for v in x) / (n - ddof)


def cusum_mean(x, confidence=0.95):
    """Standardized CUSUM of deviations from the sample mean.

    Returns ``(cusum, boundary)`` where ``cusum[k]`` is the cumulative sum of
    demeaned observations through index ``k`` divided by ``sigma * sqrt(n)``, and
    ``boundary`` is the confidence threshold. Under a stable mean the standardized
    path is a Brownian bridge (it starts and ends at 0), so the supremum of its
    absolute value is compared to the Kolmogorov critical values -- a ``cusum``
    magnitude exceeding ``boundary`` flags a structural break in the mean.
    """
    n = len(x)
    if n < 3:
        raise ValueError("need at least 3 observations")
    m = _mean(x)
    sd = math.sqrt(_var(x, ddof=1))
    if sd == 0.0:
        raise ValueError("series has zero variance")
    scale = sd * math.sqrt(n)
    cusum = []
    running = 0.0
    for v in x:
        running += (v - m)
        cusum.append(running / scale)
    # Kolmogorov critical values for sup|Brownian bridge| (a constant band).
    crit = {0.90: 1.224, 0.95: 1.358, 0.99: 1.628}.get(confidence)
    if crit is None:
        raise ValueError("confidence must be one of 0.90, 0.95, 0.99")
    return cusum, crit


def cusum_break_detected(x, confidence=0.95):
    """True if the standardized CUSUM path breaches its confidence band."""
    cusum, boundary = cusum_mean(x, confidence)
    return any(abs(c) > boundary for c in cusum)


def chow_test(x, break_index):
    """Chow F test for a break in the mean at ``break_index``.

    Compares the pooled residual sum of squares (one mean for the whole sample)
    against the sum from fitting separate means before and after the break. Under
    the no-break null the statistic is F(1, n - 2) distributed; a large value
    rejects. Returns ``(F_statistic, dof1, dof2)``.
    """
    n = len(x)
    if not (1 < break_index < n - 1):
        raise ValueError("break_index must leave >= 2 points on each side")

    def rss(seg):
        m = _mean(seg)
        return sum((v - m) ** 2 for v in seg)

    rss_pooled = rss(x)
    rss_split = rss(x[:break_index]) + rss(x[break_index:])
    k = 1                       # one parameter (the mean) per segment difference
    dof2 = n - 2 * k
    if dof2 <= 0 or rss_split <= 0.0:
        raise ValueError("degenerate configuration for the Chow test")
    f = ((rss_pooled - rss_split) / k) / (rss_split / dof2)
    return f, k, dof2
