"""Augmented Dickey-Fuller unit-root and Engle-Granger cointegration tests.

A stationary (mean-reverting) series is tradeable; a unit-root series is a random
walk. The Dickey-Fuller test regresses the change on the lagged level,

    delta y_t = alpha + rho * y_{t-1} + sum gamma_k delta y_{t-k} + eps,

and tests ``rho = 0`` (a unit root) against ``rho < 0`` (mean reversion) with the
``t`` statistic on ``rho``. Its distribution is non-standard, so it is compared to
Dickey-Fuller critical values, not the normal.

Engle-Granger cointegration then asks whether a *linear combination* of two
non-stationary series is stationary: regress ``y`` on ``x``, and ADF-test the
residual spread. If the spread has no unit root the pair is cointegrated -- the
statistical basis for a pairs trade, and the input to
:func:`~quantforge.fit_ornstein_uhlenbeck`. Pure standard library.
"""

import math

# Dickey-Fuller t critical values (constant, no trend), large-sample
# (MacKinnon). Interpolated across sample size is overkill here; these asymptotic
# values are used with a small finite-sample note.
_DF_CRIT = {0.01: -3.43, 0.05: -2.86, 0.10: -2.57}


def _ols(X, y):
    """Least squares for design matrix ``X`` (rows) and target ``y``.

    Returns (coefficients, residuals, standard_errors) via the normal equations
    solved with a Gauss-Jordan inverse.
    """
    n = len(y)
    p = len(X[0])
    # X'X and X'y.
    xtx = [[sum(X[t][i] * X[t][j] for t in range(n)) for j in range(p)]
           for i in range(p)]
    xty = [sum(X[t][i] * y[t] for t in range(n)) for i in range(p)]
    inv = _invert(xtx)
    beta = [sum(inv[i][j] * xty[j] for j in range(p)) for i in range(p)]
    resid = [y[t] - sum(X[t][i] * beta[i] for i in range(p)) for t in range(n)]
    dof = n - p
    if dof <= 0:
        raise ValueError("not enough observations for the regression")
    s2 = sum(r * r for r in resid) / dof
    se = [math.sqrt(s2 * inv[i][i]) if inv[i][i] > 0 else float("inf")
          for i in range(p)]
    return beta, resid, se


def _invert(a):
    n = len(a)
    m = [list(map(float, row)) + [1.0 if i == j else 0.0 for j in range(n)]
         for i, row in enumerate(a)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-14:
            raise ValueError("singular matrix in regression")
        m[col], m[piv] = m[piv], m[col]
        pv = m[col][col]
        m[col] = [v / pv for v in m[col]]
        for r in range(n):
            if r != col:
                f = m[r][col]
                m[r] = [m[r][k] - f * m[col][k] for k in range(2 * n)]
    return [row[n:] for row in m]


def adf_test(y, lags=0):
    """Augmented Dickey-Fuller test statistic for a unit root (constant, no trend).

    Parameters
    ----------
    y : sequence of float
        The series to test.
    lags : int
        Number of lagged differences to include (the ``augmented`` part), to soak
        up serial correlation in the residuals.

    Returns
    -------
    dict
        ``{"statistic", "rho", "reject_1pct", "reject_5pct", "reject_10pct"}``.
        A more negative statistic is stronger evidence of stationarity; the
        booleans compare it to the Dickey-Fuller critical values.
    """
    n = len(y)
    if lags < 0:
        raise ValueError("lags must be non-negative")
    if n < lags + 4:
        raise ValueError("series too short for the requested lags")

    dy = [y[i] - y[i - 1] for i in range(1, n)]
    # Build the regression delta y_t = alpha + rho y_{t-1} + sum gamma_k dy_{t-k}.
    start = lags
    X = []
    target = []
    for t in range(start, len(dy)):
        row = [1.0, y[t]]                       # intercept, lagged level y_{t-1}
        for k in range(1, lags + 1):
            row.append(dy[t - k])
        X.append(row)
        target.append(dy[t])
    beta, resid, se = _ols(X, target)
    rho = beta[1]
    stat = rho / se[1]
    return {
        "statistic": stat,
        "rho": rho,
        "reject_1pct": stat < _DF_CRIT[0.01],
        "reject_5pct": stat < _DF_CRIT[0.05],
        "reject_10pct": stat < _DF_CRIT[0.10],
    }


def engle_granger(y, x, lags=0):
    """Engle-Granger cointegration test between two series.

    Regresses ``y`` on ``x`` (with an intercept), then runs :func:`adf_test` on
    the residual spread. Rejecting the unit root in the residual means ``y`` and
    ``x`` are cointegrated.

    Returns
    -------
    dict
        ``{"hedge_ratio", "intercept", "adf", "cointegrated_5pct"}`` where ``adf``
        is the residual ADF result and ``cointegrated_5pct`` is its 5% rejection.
    """
    if len(y) != len(x):
        raise ValueError("y and x must have the same length")
    n = len(y)
    if n < 4:
        raise ValueError("need at least 4 observations")

    X = [[1.0, x[i]] for i in range(n)]
    beta, resid, _ = _ols(X, list(y))
    adf = adf_test(resid, lags=lags)
    return {
        "intercept": beta[0],
        "hedge_ratio": beta[1],
        "adf": adf,
        "cointegrated_5pct": adf["reject_5pct"],
    }
