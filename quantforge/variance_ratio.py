"""Lo-MacKinlay variance-ratio test for a random walk.

Under a random walk the variance of returns scales linearly with the horizon, so
the variance of the ``q``-period return is ``q`` times the variance of the
one-period return. The variance ratio

    VR(q) = Var(r_t(q)) / (q * Var(r_t(1)))

is therefore 1 under the random-walk null. Departures are informative:

  * ``VR(q) > 1`` -- positive serial correlation (trending / momentum);
  * ``VR(q) < 1`` -- negative serial correlation (mean reversion).

Lo and MacKinlay (1988) give a test statistic. Their heteroskedasticity-robust
``z`` divides ``VR(q) - 1`` by a standard error built from the autocovariances of
squared returns, so it stays valid under conditional heteroskedasticity (e.g.
GARCH). Under the null ``z`` is asymptotically standard normal. Pure standard
library.
"""

import math


def _mean(x):
    return sum(x) / len(x)


def variance_ratio(returns, q):
    """Lo-MacKinlay variance ratio ``VR(q)`` of a return series.

    Uses overlapping ``q``-period returns and the unbiased scaling factors from
    Lo-MacKinlay (1988). Returns 1 under a random walk, >1 for trending series,
    <1 for mean-reverting series.
    """
    n = len(returns)
    if q < 2:
        raise ValueError("q must be >= 2")
    if n <= q:
        raise ValueError("need more returns than q")

    mu = _mean(returns)
    # One-period variance (unbiased).
    var1 = sum((r - mu) ** 2 for r in returns) / (n - 1)
    if var1 == 0.0:
        raise ValueError("returns have zero variance")

    # Overlapping q-period returns, with the Lo-MacKinlay bias correction.
    m = q * (n - q + 1) * (1.0 - q / n)
    acc = 0.0
    for t in range(q - 1, n):
        s = sum(returns[t - i] for i in range(q))
        acc += (s - q * mu) ** 2
    varq = acc / m
    return varq / var1


def variance_ratio_zstat(returns, q):
    """Heteroskedasticity-robust Lo-MacKinlay ``z`` statistic for ``VR(q) = 1``.

    Divides ``VR(q) - 1`` by the robust standard error assembled from the
    autocorrelations of squared demeaned returns. Asymptotically standard normal
    under the random-walk null; ``|z| > 1.96`` rejects at 5%.
    """
    n = len(returns)
    if q < 2:
        raise ValueError("q must be >= 2")
    if n <= q:
        raise ValueError("need more returns than q")

    mu = _mean(returns)
    dev = [r - mu for r in returns]
    var1 = sum(d * d for d in dev) / n
    if var1 == 0.0:
        raise ValueError("returns have zero variance")

    vr = variance_ratio(returns, q)

    # Heteroskedasticity-robust variance of VR (Lo-MacKinlay 1988, eq. for M2).
    theta = 0.0
    for j in range(1, q):
        # delta_j: normalized autocovariance of squared deviations at lag j.
        num = 0.0
        for t in range(j, n):
            num += (dev[t] ** 2) * (dev[t - j] ** 2)
        denom = (sum(d * d for d in dev)) ** 2 / n
        delta_j = num / denom
        weight = 2.0 * (q - j) / q
        theta += weight * weight * delta_j

    se = math.sqrt(theta / n)
    if se == 0.0:
        raise ValueError("degenerate standard error")
    return (vr - 1.0) / se
