"""One-dimensional Gaussian mixture model fit by expectation-maximization.

Models a sample as a mixture of ``k`` normal components,
``p(x) = sum_j w_j N(x; mu_j, sigma_j^2)``, and fits the weights, means and variances
by the EM algorithm: the E-step computes each point's responsibility (posterior
membership) to every component, the M-step updates the parameters as
responsibility-weighted moments. The log-likelihood increases monotonically to a
local optimum. Useful for return-regime clustering (e.g. calm vs turbulent). Pure
standard library.
"""

import math


def _normal_pdf(x, mu, var):
    return math.exp(-0.5 * (x - mu) ** 2 / var) / math.sqrt(2.0 * math.pi * var)


def fit_gaussian_mixture(data, k=2, max_iter=200, tol=1e-8, seed=1234567):
    """Fit a ``k``-component 1-D Gaussian mixture by EM.

    Returns a dict with ``weights`` (mixing proportions summing to one), ``means``,
    ``variances``, ``log_likelihood`` (of the final fit) and ``n_iter``. Components
    are initialized at spread-out data quantiles (deterministic given ``seed``, which
    only jitters the initial means). The log-likelihood is non-decreasing across
    iterations. Requires at least ``k`` distinct points.
    """
    n = len(data)
    if n < k:
        raise ValueError("need at least k data points")
    if k < 1:
        raise ValueError("k must be at least 1")

    lo, hi = min(data), max(data)
    if hi <= lo:
        raise ValueError("data has zero spread")
    mean = sum(data) / n
    var0 = sum((x - mean) ** 2 for x in data) / n

    # Initialize means at evenly spaced quantiles, with a tiny seed-based jitter.
    state = seed & 0x7FFFFFFF
    def rnd():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return (state + 0.5) / 0x80000000 - 0.5
    means = [lo + (hi - lo) * (j + 1) / (k + 1) + 0.01 * (hi - lo) * rnd()
             for j in range(k)]
    variances = [var0] * k
    weights = [1.0 / k] * k

    prev_ll = -float("inf")
    n_iter = 0
    for it in range(max_iter):
        n_iter = it + 1
        # E-step: responsibilities r[i][j].
        ll = 0.0
        resp = []
        for x in data:
            comp = [weights[j] * _normal_pdf(x, means[j], variances[j])
                    for j in range(k)]
            tot = sum(comp)
            if tot <= 0.0:
                tot = 1e-300
            resp.append([c / tot for c in comp])
            ll += math.log(tot)
        # M-step.
        for j in range(k):
            nj = sum(resp[i][j] for i in range(n))
            if nj <= 1e-300:
                continue
            weights[j] = nj / n
            mj = sum(resp[i][j] * data[i] for i in range(n)) / nj
            means[j] = mj
            vj = sum(resp[i][j] * (data[i] - mj) ** 2 for i in range(n)) / nj
            variances[j] = max(vj, 1e-12)      # floor to avoid collapse
        if ll - prev_ll < tol:
            break
        prev_ll = ll

    return {"weights": weights, "means": means, "variances": variances,
            "log_likelihood": ll, "n_iter": n_iter}
