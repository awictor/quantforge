"""Mean-variance portfolio optimization from a covariance matrix.

Closed-form long/short weights for the classic problems, plus an iterative
risk-parity solver. Covariances are plain nested lists (``cov[i][j]``); no
NumPy. The analytic solutions use a small Gauss-Jordan inverse.

  * ``min_variance_weights``  -- the global minimum-variance portfolio
    ``w = C^{-1} 1 / (1^T C^{-1} 1)``;
  * ``max_sharpe_weights``    -- the tangency portfolio maximizing the Sharpe
    ratio, ``w propto C^{-1} (mu - rf)`` normalized to sum 1;
  * ``risk_parity_weights``   -- weights equalizing each asset's risk
    contribution, by fixed-point iteration;
  * ``portfolio_variance`` / ``portfolio_return`` helpers.
"""

import math
from typing import Sequence


def _matvec(A, x):
    return [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]


def _invert(matrix):
    """Gauss-Jordan inverse of a square matrix (list of lists)."""
    n = len(matrix)
    a = [list(map(float, row)) + [1.0 if i == j else 0.0 for j in range(n)]
         for i, row in enumerate(matrix)]
    for col in range(n):
        # Partial pivot.
        piv = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[piv][col]) < 1e-14:
            raise ValueError("covariance matrix is singular")
        a[col], a[piv] = a[piv], a[col]
        pv = a[col][col]
        a[col] = [v / pv for v in a[col]]
        for r in range(n):
            if r != col:
                factor = a[r][col]
                a[r] = [a[r][k] - factor * a[col][k] for k in range(2 * n)]
    return [row[n:] for row in a]


def _check_cov(cov):
    n = len(cov)
    if n == 0 or any(len(row) != n for row in cov):
        raise ValueError("cov must be a non-empty square matrix")
    return n


def portfolio_variance(weights, cov) -> float:
    """Portfolio variance ``w^T C w``."""
    cw = _matvec(cov, list(weights))
    return sum(w * v for w, v in zip(weights, cw))


def portfolio_return(weights, mean_returns) -> float:
    """Expected portfolio return ``w^T mu``."""
    return sum(w * m for w, m in zip(weights, mean_returns))


def min_variance_weights(cov) -> list:
    """Global minimum-variance weights ``C^{-1} 1 / (1^T C^{-1} 1)``.

    Fully invested (weights sum to 1); may be long/short. Requires a
    non-singular covariance matrix.
    """
    n = _check_cov(cov)
    inv = _invert(cov)
    ones = [1.0] * n
    z = _matvec(inv, ones)
    total = sum(z)
    if abs(total) < 1e-14:
        raise ValueError("degenerate covariance (zero total inverse weight)")
    return [zi / total for zi in z]


def max_sharpe_weights(mean_returns, cov, risk_free=0.0) -> list:
    """Tangency (max-Sharpe) weights ``C^{-1} (mu - rf) / sum(...)``.

    Maximizes the portfolio Sharpe ratio over fully-invested long/short
    weights. Requires the excess returns not to be orthogonal to ``C^{-1} 1``.
    """
    n = _check_cov(cov)
    if len(mean_returns) != n:
        raise ValueError("mean_returns length must match cov")
    inv = _invert(cov)
    excess = [m - risk_free for m in mean_returns]
    z = _matvec(inv, excess)
    total = sum(z)
    if abs(total) < 1e-14:
        raise ValueError("excess returns give a zero-sum tangency solution")
    return [zi / total for zi in z]


def target_return_weights(mean_returns, cov, target) -> list:
    """Minimum-variance weights achieving an exact expected return ``target``.

    Solves ``min w^T C w`` subject to ``w^T 1 = 1`` and ``w^T mu = target`` by
    the two-constraint Lagrangian. With the efficient-frontier scalars
    ``A = 1^T C^{-1} 1``, ``B = 1^T C^{-1} mu``, ``C2 = mu^T C^{-1} mu`` and
    ``D = A C2 - B^2``, the weights are

        w = C^{-1} [ (C2 - B target)/D * 1 + (A target - B)/D * mu ].

    Fully invested; may be long/short. Sweeping ``target`` traces the efficient
    frontier.
    """
    n = _check_cov(cov)
    if len(mean_returns) != n:
        raise ValueError("mean_returns length must match cov")
    inv = _invert(cov)
    ones = [1.0] * n
    mu = list(mean_returns)
    inv1 = _matvec(inv, ones)
    invmu = _matvec(inv, mu)
    A = sum(inv1)
    B = sum(mu[i] * inv1[i] for i in range(n))
    C2 = sum(mu[i] * invmu[i] for i in range(n))
    D = A * C2 - B * B
    if abs(D) < 1e-14:
        raise ValueError("degenerate frontier (returns collinear with ones)")
    g = (C2 - B * target) / D
    h = (A * target - B) / D
    return [g * inv1[i] + h * invmu[i] for i in range(n)]


def efficient_frontier(mean_returns, cov, targets) -> list:
    """Efficient frontier as ``(target_return, portfolio_std)`` pairs.

    For each requested expected return in ``targets`` solves
    :func:`target_return_weights` and reports the achieved return with the
    portfolio standard deviation ``sqrt(w^T C w)``.
    """
    out = []
    for tgt in targets:
        w = target_return_weights(mean_returns, cov, tgt)
        out.append((tgt, math.sqrt(portfolio_variance(w, cov))))
    return out


def risk_parity_weights(cov, tol=1e-10, max_iter=1000) -> list:
    """Equal-risk-contribution (risk-parity) weights.

    Solves for positive weights whose marginal risk contributions
    ``w_i (C w)_i`` are equal, by the standard fixed-point iteration
    ``w_i <- (target / (C w)_i)`` renormalized -- convergent for a positive-
    definite covariance. Weights are long-only and sum to 1.
    """
    n = _check_cov(cov)
    w = [1.0 / n] * n
    for _ in range(max_iter):
        cw = _matvec(cov, w)
        if any(v <= 0.0 for v in cw):
            raise ValueError("covariance not positive-definite for risk parity")
        # ERC fixed point with sqrt damping: w_i <- w_i sqrt(target / RC_i),
        # RC_i = w_i (C w)_i, so w_i sqrt(1/(w_i cw_i)) = sqrt(w_i / cw_i). The
        # damping avoids the oscillation the raw 1/(Cw) update shows on diagonal
        # covariances; the fixed point equalizes risk contributions (w_i ∝
        # 1/sigma_i for a diagonal cov).
        new = [math.sqrt(w[i] / cw[i]) for i in range(n)]
        ssum = sum(new)
        new = [v / ssum for v in new]
        if max(abs(new[i] - w[i]) for i in range(n)) < tol:
            w = new
            break
        w = new
    return w
