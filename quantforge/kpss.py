"""KPSS test for stationarity (Kwiatkowski-Phillips-Schmidt-Shin).

Where the augmented Dickey-Fuller test takes a unit root as the null, KPSS takes
*stationarity* as the null -- so the two are complementary confirmatory tests. The
statistic is

    eta = (1 / n^2) sum_t S_t^2 / s^2(l),

where ``S_t`` is the cumulative sum of the residuals from regressing the series on a
constant (level stationarity) or a constant plus trend (trend stationarity), and
``s^2(l)`` is the Newey-West long-run variance of the residuals. A large statistic
rejects stationarity (evidence of a unit root / drift). Compared against the
Kwiatkowski et al. (1992) asymptotic critical values. Pure standard library.
"""

from .hac import newey_west_variance

# Asymptotic critical values (Kwiatkowski et al. 1992, Table 1).
_CRIT_LEVEL = {0.10: 0.347, 0.05: 0.463, 0.025: 0.574, 0.01: 0.739}
_CRIT_TREND = {0.10: 0.119, 0.05: 0.146, 0.025: 0.176, 0.01: 0.216}


def kpss_test(y, lags=None, regression="c"):
    """KPSS stationarity test statistic and approximate p-value.

    ``regression`` is ``"c"`` for level stationarity (residuals about the mean) or
    ``"ct"`` for trend stationarity (residuals about a fitted linear trend).
    ``lags`` sets the Newey-West bandwidth for the long-run variance (defaults to
    ``floor(4 (n/100)^{1/4})``). Returns ``(eta, p_value)``; a *small* p-value
    rejects stationarity (unlike ADF, where a small p-value supports it). The p-value
    is interpolated/clamped against the asymptotic critical values, so it is reported
    within ``[0.01, 0.10]`` at the bounds.
    """
    n = len(y)
    if n < 4:
        raise ValueError("need at least 4 observations")
    if regression not in ("c", "ct"):
        raise ValueError("regression must be 'c' or 'ct'")
    if lags is None:
        lags = int(4.0 * (n / 100.0) ** 0.25)

    # Residuals from the deterministic regression.
    if regression == "c":
        mean = sum(y) / n
        resid = [v - mean for v in y]
    else:
        xs = list(range(n))
        mx = (n - 1) / 2.0
        my = sum(y) / n
        sxx = sum((x - mx) ** 2 for x in xs)
        sxy = sum((xs[i] - mx) * (y[i] - my) for i in range(n))
        slope = sxy / sxx
        intercept = my - slope * mx
        resid = [y[i] - (slope * i + intercept) for i in range(n)]

    # Cumulative sum of residuals.
    cumulative = []
    acc = 0.0
    for r in resid:
        acc += r
        cumulative.append(acc)
    s_sq = sum(c * c for c in cumulative) / (n * n)
    lr_var = newey_west_variance(resid, lags)
    if lr_var <= 0.0:
        raise ValueError("non-positive long-run variance")
    eta = s_sq / lr_var

    crit = _CRIT_LEVEL if regression == "c" else _CRIT_TREND
    p = _kpss_pvalue(eta, crit)
    return eta, p


def _kpss_pvalue(eta, crit):
    """Interpolate a p-value from the KPSS critical-value table (clamped)."""
    levels = sorted(crit)                       # ascending p: 0.01, 0.025, 0.05, 0.10
    values = [crit[p] for p in levels]          # descending critical values
    # Larger eta -> smaller p. Clamp outside the table.
    if eta >= values[0]:                        # >= 1% critical value
        return 0.01
    if eta <= values[-1]:                       # <= 10% critical value
        return 0.10
    # Linear interpolation between adjacent (value, level) points.
    for i in range(len(values) - 1):
        hi_v, lo_v = values[i], values[i + 1]
        if lo_v <= eta <= hi_v:
            hi_p, lo_p = levels[i], levels[i + 1]
            w = (eta - lo_v) / (hi_v - lo_v)
            return lo_p + w * (hi_p - lo_p)
    return 0.10
