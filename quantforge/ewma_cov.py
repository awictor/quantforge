"""EWMA (RiskMetrics) covariance and correlation matrices.

The exponentially-weighted moving-average covariance weights recent
co-movements more heavily, tracking a time-varying covariance far better than an
equal-weighted sample. For a decay ``lambda`` in ``(0, 1)`` the recursion is

    Sigma_t = lambda Sigma_{t-1} + (1 - lambda) r_t r_t',

seeded from the sample covariance. RiskMetrics uses ``lambda = 0.94`` for daily
data. The estimate is symmetric and positive semidefinite. Pure standard library.
"""


def _sample_cov(returns):
    n = len(returns)
    p = len(returns[0])
    mean = [sum(returns[t][i] for t in range(n)) / n for i in range(p)]
    cov = [[0.0] * p for _ in range(p)]
    for t in range(n):
        dev = [returns[t][i] - mean[i] for i in range(p)]
        for i in range(p):
            for j in range(p):
                cov[i][j] += dev[i] * dev[j]
    for i in range(p):
        for j in range(p):
            cov[i][j] /= n
    return cov


def ewma_covariance_matrix(returns, lam=0.94):
    """EWMA covariance matrix of a return panel (rows = periods, cols = assets).

    Recurses ``Sigma_t = lam Sigma_{t-1} + (1-lam) r_t r_t'`` from the sample
    covariance seed, treating returns as mean-zero (the RiskMetrics convention).
    Returns the final ``p x p`` symmetric positive-semidefinite matrix.
    """
    n = len(returns)
    if n < 2:
        raise ValueError("need at least 2 observations")
    p = len(returns[0])
    if any(len(r) != p for r in returns):
        raise ValueError("all rows must have the same length")
    if not (0.0 < lam < 1.0):
        raise ValueError("lam must be in (0, 1)")

    sigma = _sample_cov(returns)
    for t in range(n):
        r = returns[t]
        for i in range(p):
            for j in range(p):
                sigma[i][j] = lam * sigma[i][j] + (1.0 - lam) * r[i] * r[j]
    # Symmetrize against floating error.
    for i in range(p):
        for j in range(i + 1, p):
            m = 0.5 * (sigma[i][j] + sigma[j][i])
            sigma[i][j] = sigma[j][i] = m
    return sigma


def ewma_correlation_matrix(returns, lam=0.94):
    """EWMA correlation matrix: the EWMA covariance normalized by its diagonal.

    Unit diagonal, off-diagonals in ``[-1, 1]``.
    """
    cov = ewma_covariance_matrix(returns, lam)
    p = len(cov)
    sd = [cov[i][i] ** 0.5 for i in range(p)]
    corr = [[0.0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            denom = sd[i] * sd[j]
            c = cov[i][j] / denom if denom > 0 else (1.0 if i == j else 0.0)
            corr[i][j] = max(-1.0, min(1.0, c))     # clamp floating error
    return corr
