"""Mack's (1993) distribution-free standard error for chain-ladder reserves.

The chain-ladder point estimate says nothing about how uncertain the reserve is.
Mack (1993, ASTIN Bulletin 23(2)) derived the mean-squared error of prediction of
the reserve under three assumptions -- linearity of expected development, uncorrelated
development factors, and a variance proportional to the prior cumulative claims,

    E[C_{i,k+1} | C_{i,k}] = f_k C_{i,k},
    Var[C_{i,k+1} | C_{i,k}] = sigma_k^2 C_{i,k}.

The age-to-age factors ``f_k`` are the usual volume-weighted estimates; the variance
parameters are

    sigma_k^2 = (1 / (n-k-2)) sum_i C_{i,k} (C_{i,k+1}/C_{i,k} - f_k)^2.

The MSEP of a single accident-year ultimate combines a process-error and an
estimation-error term,

    mse(hat C_{i,I}) = hat C_{i,I}^2 sum_{k=I-i}^{I-1} (sigma_k^2 / f_k^2)
                       (1 / hat C_{i,k} + 1 / S_k),   S_k = sum_j C_{j,k},

and the total-reserve MSEP adds the correlation between accident years that share
the same estimated factors. The last variance parameter, which cannot be estimated
from a single residual, is extrapolated by Mack's log-linear rule. Pure standard
library; validated against the Taylor-Ashe triangle in Mack's paper (total standard
error 2,447,095 on a reserve of 18,680,856).
"""


def _dev_factors_and_S(triangle):
    """Volume-weighted factors ``f_k`` and column sums ``S_k = sum_i C_{i,k}``.

    ``S_k`` runs over the accident years that contribute to ``f_k`` (those observed
    at both ages ``k`` and ``k+1``), matching the denominator of the factor.
    """
    n = len(triangle)
    f = []
    S = []
    for k in range(n - 1):
        num = 0.0
        den = 0.0
        for i in range(n):
            if len(triangle[i]) > k + 1:
                den += triangle[i][k]
                num += triangle[i][k + 1]
        if den <= 0.0:
            raise ValueError("no data to estimate a development factor")
        f.append(num / den)
        S.append(den)
    return f, S


def _sigma2(triangle, f):
    """Mack variance parameters ``sigma_k^2``, last one extrapolated.

    ``sigma_k^2 = (1/(n-k-2)) sum_i C_{i,k} (C_{i,k+1}/C_{i,k} - f_k)^2`` for
    ``k <= n-3``. The final parameter ``sigma_{n-2}^2`` has only one observation, so
    Mack's rule sets it to ``min(sigma_{n-3}^4 / sigma_{n-4}^2, sigma_{n-3}^2,
    sigma_{n-4}^2)`` -- a log-linear extrapolation clamped to the last two values.
    """
    n = len(triangle)
    sigma2 = []
    for k in range(n - 1):
        m = 0  # number of accident years observed at both k and k+1
        acc = 0.0
        for i in range(n):
            if len(triangle[i]) > k + 1:
                m += 1
                ratio = triangle[i][k + 1] / triangle[i][k]
                acc += triangle[i][k] * (ratio - f[k]) ** 2
        if m >= 2:
            sigma2.append(acc / (m - 1))
        else:
            # Only one residual (k = n-2): extrapolate from the previous two.
            if len(sigma2) >= 2 and sigma2[-2] > 0.0:
                a, b = sigma2[-1], sigma2[-2]  # sigma_{n-3}^2, sigma_{n-4}^2
                sigma2.append(min(a * a / b, a, b))
            elif sigma2:
                sigma2.append(sigma2[-1])
            else:
                sigma2.append(0.0)
    return sigma2


def mack_standard_error(triangle):
    """Chain-ladder reserves with Mack's prediction standard errors.

    ``triangle`` is a cumulative run-off triangle: ``triangle[i]`` holds the observed
    cumulative claims for accident year ``i`` at development ages ``0 ..``, with the
    most recent accident year the shortest row. Returns a dict with

    - ``factors``, ``sigma2``: estimated age-to-age factors and variance parameters,
    - ``ultimate``, ``reserve``: per accident year,
    - ``std_error``: Mack standard error of each accident-year reserve,
    - ``cv``: coefficient of variation ``std_error / reserve`` (0 for a zero reserve),
    - ``total_reserve``, ``total_std_error``, ``total_cv``: for the reserve sum,
      including the between-year correlation term.

    Raises ``ValueError`` for fewer than two accident years or non-positive entries.
    """
    n = len(triangle)
    if n < 2:
        raise ValueError("need at least 2 accident years")
    for row in triangle:
        if any(c <= 0.0 for c in row):
            raise ValueError("cumulative claims must be positive")

    f, S = _dev_factors_and_S(triangle)
    sigma2 = _sigma2(triangle, f)

    # Project each row to ultimate, keeping the full developed path C_{i,k}.
    full = []
    for i in range(n):
        row = list(triangle[i])
        for k in range(len(row) - 1, n - 1):
            row.append(row[k] * f[k])
        full.append(row)
    ultimate = [full[i][n - 1] for i in range(n)]
    latest = [triangle[i][-1] for i in range(n)]
    reserve = [ultimate[i] - latest[i] for i in range(n)]

    # Per-year mean-squared error of prediction (Mack eq. 3).
    mse = [0.0] * n
    for i in range(n):
        start = len(triangle[i]) - 1  # first age to be developed
        s = 0.0
        for k in range(start, n - 1):
            term = (sigma2[k] / (f[k] * f[k])) * (1.0 / full[i][k] + 1.0 / S[k])
            s += term
        mse[i] = ultimate[i] * ultimate[i] * s
    std_error = [msei ** 0.5 for msei in mse]
    cv = [std_error[i] / reserve[i] if reserve[i] > 0.0 else 0.0 for i in range(n)]

    # Total reserve MSEP: sum of per-year MSEP plus the correlation between years
    # that share estimated factors (Mack eq. 7).
    total_mse = sum(mse)
    for i in range(n):
        tail_sum = sum(ultimate[j] for j in range(i + 1, n))
        if tail_sum == 0.0:
            continue
        start = len(triangle[i]) - 1
        cross = 0.0
        for k in range(start, n - 1):
            cross += (sigma2[k] / (f[k] * f[k])) / S[k]
        total_mse += 2.0 * ultimate[i] * tail_sum * cross

    total_reserve = sum(reserve)
    total_std_error = total_mse ** 0.5
    total_cv = total_std_error / total_reserve if total_reserve > 0.0 else 0.0

    return {
        "factors": f,
        "sigma2": sigma2,
        "ultimate": ultimate,
        "reserve": reserve,
        "std_error": std_error,
        "cv": cv,
        "total_reserve": total_reserve,
        "total_std_error": total_std_error,
        "total_cv": total_cv,
    }
