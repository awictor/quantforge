"""Diebold-Mariano test of equal predictive accuracy.

Compares two competing forecasts of the same target by their loss differential
``d_t = L(e1_t) - L(e2_t)`` (default squared-error loss ``L(e) = e^2``). Under the
null of equal accuracy ``E[d_t] = 0``. The test statistic standardizes the mean
loss differential by its Newey-West (HAC) standard error, which accounts for the
serial correlation that forecast errors of horizon ``h`` inevitably carry:

    DM = mean(d) / sqrt( long_run_var(d) / n ),

asymptotically ``N(0, 1)``. A significantly negative ``DM`` means forecast 1 is more
accurate (smaller loss); positive means forecast 2 wins. The Harvey-Leybourne-
Newbold small-sample correction is applied and a Student-t reference used. Pure
standard library on top of the HAC variance and the t distribution.
"""

import math

from .hac import newey_west_variance
from .student_t import t_cdf


def diebold_mariano(errors1, errors2, h=1, power=2):
    """Diebold-Mariano statistic and two-sided p-value for equal predictive accuracy.

    ``errors1`` / ``errors2`` are the forecast-error series of the two models
    (target minus forecast), ``h`` the forecast horizon (sets the number of
    autocovariance lags ``h - 1`` in the HAC variance), and ``power`` the loss
    exponent (2 = squared error, 1 = absolute error). Returns ``(DM, p_value)`` with
    the Harvey-Leybourne-Newbold small-sample correction and a Student-t reference on
    ``n - 1`` degrees of freedom. A negative ``DM`` favors the first forecast.
    """
    n = len(errors1)
    if n < 2 or len(errors2) != n:
        raise ValueError("error series must be equal-length with at least 2 points")
    if h < 1:
        raise ValueError("h must be at least 1")
    d = [abs(errors1[i]) ** power - abs(errors2[i]) ** power for i in range(n)]
    mean_d = sum(d) / n
    lags = h - 1
    lr_var = newey_west_variance(d, lags)
    if lr_var <= 0.0:
        raise ValueError("non-positive long-run variance (identical forecasts?)")
    dm = mean_d / math.sqrt(lr_var / n)
    # Harvey-Leybourne-Newbold small-sample correction.
    correction = math.sqrt((n + 1 - 2 * h + h * (h - 1) / n) / n)
    dm_corrected = dm * correction
    p = 2.0 * (1.0 - t_cdf(abs(dm_corrected), n - 1))
    return dm_corrected, p
