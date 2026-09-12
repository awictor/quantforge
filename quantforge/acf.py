"""Autocorrelation and partial-autocorrelation functions (ACF / PACF).

These are the model-identification workhorses of time-series analysis. The ACF at
lag ``k`` is the correlation between the series and its ``k``-lagged self; the PACF
is the correlation after removing the linear effect of the intervening lags. Their
signatures identify ARMA order:

  * an AR(p) process has a PACF that cuts off after lag ``p`` and an ACF that
    decays geometrically;
  * an MA(q) process has an ACF that cuts off after lag ``q``.

The PACF is computed by the Durbin-Levinson recursion, which solves the Yule-Walker
equations level by level; the last reflection coefficient at each order is that
order's partial autocorrelation. Pure standard library.
"""


def _autocov(x, k, mean):
    n = len(x)
    return sum((x[t] - mean) * (x[t - k] - mean) for t in range(k, n)) / n


def acf(x, nlags=20):
    """Autocorrelation function up to ``nlags`` (lag 0 included, always 1).

    Uses the biased (divisor ``n``) autocovariance, giving a positive-semidefinite
    sequence. Returns a list of length ``nlags + 1``.
    """
    n = len(x)
    if n < 2:
        raise ValueError("need at least 2 points")
    if nlags < 1 or nlags >= n:
        raise ValueError("nlags must satisfy 1 <= nlags < n")
    mean = sum(x) / n
    g0 = _autocov(x, 0, mean)
    if g0 == 0.0:
        raise ValueError("series has zero variance")
    return [_autocov(x, k, mean) / g0 for k in range(nlags + 1)]


def pacf(x, nlags=20):
    """Partial autocorrelation function via the Durbin-Levinson recursion.

    Returns a list of length ``nlags + 1`` with ``pacf[0] = 1`` by convention and
    ``pacf[k]`` the order-``k`` partial autocorrelation. For an AR(p) process the
    PACF is ~0 beyond lag ``p``.
    """
    n = len(x)
    if nlags < 1 or nlags >= n:
        raise ValueError("nlags must satisfy 1 <= nlags < n")
    r = acf(x, nlags)             # r[0]=1, r[1..nlags]
    phi = [[0.0] * (nlags + 1) for _ in range(nlags + 1)]
    out = [1.0] + [0.0] * nlags
    # Durbin-Levinson.
    phi[1][1] = r[1]
    out[1] = r[1]
    for k in range(2, nlags + 1):
        num = r[k] - sum(phi[k - 1][j] * r[k - j] for j in range(1, k))
        den = 1.0 - sum(phi[k - 1][j] * r[j] for j in range(1, k))
        phi_kk = num / den if den != 0.0 else 0.0
        phi[k][k] = phi_kk
        for j in range(1, k):
            phi[k][j] = phi[k - 1][j] - phi_kk * phi[k - 1][k - j]
        out[k] = phi_kk
    return out
