"""Engle's ARCH-LM test for conditional heteroskedasticity (volatility clustering).

Financial returns show volatility clustering: large moves follow large moves. Engle
(1982) tests for it by regressing the squared residuals on their own lags,

    e_t^2 = a_0 + a_1 e_{t-1}^2 + ... + a_q e_{t-q}^2 + u_t,

under the null of no ARCH effect all slopes are zero, so ``LM = n * R^2`` is
asymptotically chi-square with ``q`` degrees of freedom. A small p-value means the
variance depends on the recent past -- a GARCH-type model is warranted. Pure
standard library on top of the OLS fit and the chi-square distribution.
"""

from .ols import ols_fit
from .distributions import chi2_sf


def arch_lm_test(residuals, lags=1):
    """Engle ARCH-LM test on a residual (or return) series.

    Regresses the squared series on ``lags`` of its own past and returns
    ``(LM, p_value)`` with ``LM = m * R^2`` (``m`` the number of regression rows)
    referenced to a chi-square with ``lags`` degrees of freedom. A small p-value
    rejects "no ARCH effect", i.e. detects volatility clustering. Subtracts the mean
    first, so it works on returns directly. Requires ``len > 2 * lags + 1``.
    """
    n = len(residuals)
    if lags < 1:
        raise ValueError("lags must be at least 1")
    if n <= 2 * lags + 1:
        raise ValueError("series too short for the chosen number of lags")
    mean = sum(residuals) / n
    e2 = [(r - mean) ** 2 for r in residuals]

    # Build the lagged design: target e2[t], regressors e2[t-1..t-lags].
    y = []
    X = []
    for t in range(lags, n):
        y.append(e2[t])
        X.append([e2[t - k] for k in range(1, lags + 1)])
    m = len(y)
    fit = ols_fit(X, y)
    r_squared = fit["r_squared"]
    lm = m * r_squared
    return lm, chi2_sf(lm, lags)
