"""Scoring rules for probabilistic forecasts.

Point-forecast metrics ignore uncertainty; these score quantile, interval and full
distributional forecasts:

- ``pinball_loss`` -- the quantile (pinball) loss, minimized in expectation at the
  true ``tau``-quantile; the loss behind quantile regression,
- ``interval_score`` -- the Winkler / Gneiting-Raftery score for a central
  prediction interval, rewarding narrow intervals but penalizing misses,
- ``coverage`` -- the empirical fraction of actuals inside their intervals (should
  match the nominal level),
- ``crps_ensemble`` -- the continuous ranked probability score of an ensemble
  (sample) forecast, the distributional generalization of absolute error.

All are proper scoring rules (lower is better). Pure standard library.
"""


def pinball_loss(actual, quantile_forecast, tau):
    """Average pinball (quantile) loss at level ``tau``.

    ``L = mean( tau (a - q)      if a >= q
                (1 - tau)(q - a)  otherwise )``.
    Minimized in expectation when ``quantile_forecast`` is the true ``tau``-quantile
    of the target. Aligned series; ``tau`` in ``(0, 1)``.
    """
    n = len(actual)
    if n == 0 or len(quantile_forecast) != n:
        raise ValueError("actual and quantile_forecast must align, non-empty")
    if not (0.0 < tau < 1.0):
        raise ValueError("tau must be in (0, 1)")
    total = 0.0
    for i in range(n):
        diff = actual[i] - quantile_forecast[i]
        total += tau * diff if diff >= 0.0 else (1.0 - tau) * (-diff)
    return total / n


def interval_score(actual, lower, upper, alpha=0.1):
    """Winkler interval score for central ``1 - alpha`` prediction intervals.

    ``S = (upper - lower) + (2/alpha)(lower - a) if a < lower
                          + (2/alpha)(a - upper) if a > upper``.
    Rewards narrow intervals and penalizes actuals falling outside, scaled so the
    penalty grows as the nominal coverage tightens. Lower is better. Aligned series.
    """
    n = len(actual)
    if n == 0 or len(lower) != n or len(upper) != n:
        raise ValueError("actual, lower, upper must align, non-empty")
    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must be in (0, 1)")
    total = 0.0
    for i in range(n):
        width = upper[i] - lower[i]
        s = width
        if actual[i] < lower[i]:
            s += (2.0 / alpha) * (lower[i] - actual[i])
        elif actual[i] > upper[i]:
            s += (2.0 / alpha) * (actual[i] - upper[i])
        total += s
    return total / n


def coverage(actual, lower, upper):
    """Empirical coverage: fraction of actuals within ``[lower, upper]``.

    Should match the interval's nominal level (e.g. ~0.9 for a 90% interval).
    """
    n = len(actual)
    if n == 0 or len(lower) != n or len(upper) != n:
        raise ValueError("actual, lower, upper must align, non-empty")
    inside = sum(1 for i in range(n) if lower[i] <= actual[i] <= upper[i])
    return inside / n


def crps_ensemble(actual, ensemble):
    """Continuous ranked probability score of an ensemble forecast (single target).

    Uses the empirical-CDF form ``CRPS = mean|X - y| - 0.5 mean|X - X'|`` where ``X``,
    ``X'`` are independent ensemble members and ``y`` the realized value. Reduces to
    the absolute error for a deterministic (single-member) forecast, and is zero when
    every member equals the target. ``ensemble`` is the list of member forecasts.
    """
    m = len(ensemble)
    if m == 0:
        raise ValueError("ensemble must be non-empty")
    term1 = sum(abs(x - actual) for x in ensemble) / m
    if m == 1:
        return term1
    term2 = 0.0
    for i in range(m):
        for j in range(m):
            term2 += abs(ensemble[i] - ensemble[j])
    term2 /= (m * m)
    return term1 - 0.5 * term2
