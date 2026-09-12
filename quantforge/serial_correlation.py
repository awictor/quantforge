"""Serial-correlation tests: Ljung-Box, Box-Pierce, Durbin-Watson.

After fitting a model, its residuals should be white noise. These portmanteau
tests check whether the first ``h`` autocorrelations are jointly zero:

  * ``box_pierce``  -- Q = n * sum_{k=1}^{h} rho_k^2.
  * ``ljung_box``   -- Q = n(n+2) * sum_{k=1}^{h} rho_k^2 / (n - k), a small-sample
    refinement. Under the white-noise null Q ~ chi-square(h), so a large Q (small
    p-value) rejects "no autocorrelation".
  * ``durbin_watson`` -- d ~ 2(1 - rho_1); near 2 means no first-order
    autocorrelation, near 0 strong positive, near 4 strong negative.

Pure standard library (includes a chi-square survival function via the regularized
incomplete gamma).
"""

import math


def _autocorr(x, k):
    n = len(x)
    m = sum(x) / n
    denom = sum((v - m) ** 2 for v in x)
    if denom == 0.0:
        raise ValueError("series has zero variance")
    num = sum((x[t] - m) * (x[t - k] - m) for t in range(k, n))
    return num / denom


def _lower_gamma_reg(s, x):
    """Regularized lower incomplete gamma P(s, x) = gamma(s, x) / Gamma(s)."""
    if x < 0.0 or s <= 0.0:
        raise ValueError("invalid arguments to incomplete gamma")
    if x == 0.0:
        return 0.0
    if x < s + 1.0:
        # Series expansion.
        term = 1.0 / s
        total = term
        n = s
        for _ in range(1000):
            n += 1.0
            term *= x / n
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        return total * math.exp(-x + s * math.log(x) - math.lgamma(s))
    # Continued fraction for the upper gamma, then complement.
    tiny = 1e-300
    b = x + 1.0 - s
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 1000):
        an = -i * (i - s)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    q = math.exp(-x + s * math.log(x) - math.lgamma(s)) * h
    return 1.0 - q


def _chi2_sf(x, df):
    """Chi-square survival function P(X > x) for ``df`` degrees of freedom."""
    if x <= 0.0:
        return 1.0
    return 1.0 - _lower_gamma_reg(df / 2.0, x / 2.0)


def box_pierce(x, lags=10):
    """Box-Pierce portmanteau statistic and p-value.

    Returns ``(Q, p_value)`` with ``Q = n * sum_{k=1}^{lags} rho_k^2`` and the
    chi-square(``lags``) p-value. Small p rejects the white-noise null.
    """
    n = len(x)
    if lags < 1 or lags >= n:
        raise ValueError("lags must satisfy 1 <= lags < n")
    q = n * sum(_autocorr(x, k) ** 2 for k in range(1, lags + 1))
    return q, _chi2_sf(q, lags)


def ljung_box(x, lags=10):
    """Ljung-Box portmanteau statistic and p-value (small-sample refinement).

    Returns ``(Q, p_value)`` with ``Q = n(n+2) sum_{k=1}^{lags} rho_k^2 / (n-k)``
    and the chi-square(``lags``) p-value. Large Q / small p rejects "the first
    ``lags`` autocorrelations are jointly zero".
    """
    n = len(x)
    if lags < 1 or lags >= n:
        raise ValueError("lags must satisfy 1 <= lags < n")
    q = n * (n + 2.0) * sum(_autocorr(x, k) ** 2 / (n - k)
                            for k in range(1, lags + 1))
    return q, _chi2_sf(q, lags)


def durbin_watson(x):
    """Durbin-Watson statistic ``d = sum (x_t - x_{t-1})^2 / sum x_t^2``.

    Approximately ``2(1 - rho_1)``: near 2 = no first-order autocorrelation, near
    0 = strong positive, near 4 = strong negative. Computed on the demeaned series.
    """
    n = len(x)
    if n < 2:
        raise ValueError("need at least 2 points")
    m = sum(x) / n
    d = [v - m for v in x]
    denom = sum(v * v for v in d)
    if denom == 0.0:
        raise ValueError("series has zero variance")
    num = sum((d[t] - d[t - 1]) ** 2 for t in range(1, n))
    return num / denom
