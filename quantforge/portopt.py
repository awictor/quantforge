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


def _norm_ppf(p):
    """Inverse standard-normal CDF (Acklam's rational approximation)."""
    if not (0.0 < p < 1.0):
        raise ValueError("p must be in (0, 1)")
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


def portfolio_var(weights, cov, mean_returns=None, confidence=0.95,
                  horizon=1.0) -> float:
    """Parametric (Gaussian) Value-at-Risk of a portfolio, as a positive loss.

    ``VaR = z * sigma_p * sqrt(horizon) - mu_p * horizon`` where ``sigma_p`` is
    the portfolio standard deviation, ``mu_p`` the expected return (0 if
    ``mean_returns`` is omitted), and ``z`` the standard-normal quantile at
    ``confidence``. Returned as a non-negative loss figure.
    """
    sd = math.sqrt(portfolio_variance(weights, cov))
    mu = 0.0 if mean_returns is None else portfolio_return(weights, mean_returns)
    z = _norm_ppf(confidence)
    return z * sd * math.sqrt(horizon) - mu * horizon


def portfolio_cvar(weights, cov, mean_returns=None, confidence=0.95,
                   horizon=1.0) -> float:
    """Parametric (Gaussian) Conditional VaR / expected shortfall, as a loss.

    ``CVaR = phi(z)/(1-c) * sigma_p * sqrt(horizon) - mu_p * horizon``, the mean
    loss beyond the VaR under normality. Exceeds :func:`portfolio_var`.
    """
    sd = math.sqrt(portfolio_variance(weights, cov))
    mu = 0.0 if mean_returns is None else portfolio_return(weights, mean_returns)
    z = _norm_ppf(confidence)
    phi = math.exp(-0.5 * z * z) / math.sqrt(2.0 * math.pi)
    return phi / (1.0 - confidence) * sd * math.sqrt(horizon) - mu * horizon


def component_var(weights, cov) -> list:
    """Component (risk-contribution) VaR: each asset's share of portfolio vol.

    The marginal contribution ``(C w)_i / sigma_p`` times ``w_i`` gives the
    component ``w_i (C w)_i / sigma_p``; the components sum to the portfolio
    standard deviation. Scale by the VaR z-quantile to get VaR contributions.
    """
    n = _check_cov(cov)
    cw = _matvec(cov, list(weights))
    sd = math.sqrt(portfolio_variance(weights, cov))
    if sd <= 0.0:
        raise ValueError("portfolio variance must be positive")
    return [weights[i] * cw[i] / sd for i in range(n)]


def implied_equilibrium_returns(cov, market_weights, risk_aversion=2.5) -> list:
    """Reverse-optimized (implied) equilibrium excess returns ``Pi = lambda C w``.

    Given the market-cap weights and a risk-aversion ``lambda``, the returns
    that make those weights mean-variance optimal are ``lambda C w`` -- the
    Black-Litterman market prior.
    """
    n = _check_cov(cov)
    if len(market_weights) != n:
        raise ValueError("market_weights length must match cov")
    cw = _matvec(cov, list(market_weights))
    return [risk_aversion * v for v in cw]


def black_litterman_returns(cov, market_weights, P, Q, tau=0.05,
                            risk_aversion=2.5, omega=None) -> list:
    """Black-Litterman posterior expected returns blending prior and views.

    The market-equilibrium prior ``Pi = lambda C w`` is combined with ``k``
    linear views ``P mu = Q`` (each row of ``P`` a portfolio, ``Q`` its expected
    return) of uncertainty ``omega`` (defaults to ``diag(tau P C P^T)``). The
    posterior mean is the standard closed form

        mu = [ (tau C)^{-1} + P^T Omega^{-1} P ]^{-1}
             [ (tau C)^{-1} Pi + P^T Omega^{-1} Q ].

    With no views (empty ``P``) it returns the prior ``Pi``.
    """
    n = _check_cov(cov)
    pi = implied_equilibrium_returns(cov, market_weights, risk_aversion)
    if not P:
        return pi
    k = len(P)
    if any(len(row) != n for row in P):
        raise ValueError("each view in P must have length n")
    if len(Q) != k:
        raise ValueError("Q length must match the number of views")
    tau_cov = [[tau * cov[i][j] for j in range(n)] for i in range(n)]
    inv_tau_cov = _invert(tau_cov)
    # Omega: default diag(tau P C P^T).
    if omega is None:
        omega_diag = []
        for r in range(k):
            cp = _matvec(cov, P[r])
            omega_diag.append(tau * sum(P[r][j] * cp[j] for j in range(n)))
        omega_inv = [[0.0] * k for _ in range(k)]
        for r in range(k):
            if omega_diag[r] <= 0.0:
                raise ValueError("view has zero prior variance")
            omega_inv[r][r] = 1.0 / omega_diag[r]
    else:
        omega_inv = _invert(omega)
    # A = inv_tau_cov + P^T Omega^{-1} P ; b = inv_tau_cov Pi + P^T Omega^{-1} Q.
    # Precompute P^T Omega^{-1}.
    pt_oi = [[sum(P[r][i] * omega_inv[r][c] for r in range(k)) for c in range(k)]
             for i in range(n)]
    ptoip = [[sum(pt_oi[i][c] * P[c][j] for c in range(k)) for j in range(n)]
             for i in range(n)]
    A = [[inv_tau_cov[i][j] + ptoip[i][j] for j in range(n)] for i in range(n)]
    itc_pi = _matvec(inv_tau_cov, pi)
    ptoi_q = [sum(pt_oi[i][c] * Q[c] for c in range(k)) for i in range(n)]
    b = [itc_pi[i] + ptoi_q[i] for i in range(n)]
    return _matvec(_invert(A), b)


def marginal_var(weights, cov, confidence=0.95, horizon=1.0) -> list:
    """Marginal VaR: sensitivity of the portfolio VaR to each weight.

    ``dVaR/dw_i = z sqrt(horizon) (C w)_i / sigma_p`` (the zero-mean parametric
    VaR). Multiplying by ``w_i`` gives the component VaR, and the dot product
    ``sum_i w_i * marginal_i`` recovers the total VaR (VaR is homogeneous of
    degree 1 in the weights).
    """
    n = _check_cov(cov)
    cw = _matvec(cov, list(weights))
    sd = math.sqrt(portfolio_variance(weights, cov))
    if sd <= 0.0:
        raise ValueError("portfolio variance must be positive")
    z = _norm_ppf(confidence) * math.sqrt(horizon)
    return [z * cw[i] / sd for i in range(n)]


def var_budget(weights, cov) -> list:
    """Percentage VaR budget: each asset's fractional share of portfolio risk.

    ``w_i (C w)_i / (w^T C w)`` -- the component VaRs normalized to sum to 1.
    Independent of the confidence level and horizon (they cancel). Shows how the
    total risk is distributed across positions; equal entries mean risk parity.
    """
    n = _check_cov(cov)
    cw = _matvec(cov, list(weights))
    total = portfolio_variance(weights, cov)
    if total <= 0.0:
        raise ValueError("portfolio variance must be positive")
    return [weights[i] * cw[i] / total for i in range(n)]


def diversification_ratio(weights, cov) -> float:
    """Diversification ratio ``(sum_i w_i sigma_i) / sqrt(w^T C w)``.

    The weighted average of the assets' standalone volatilities over the
    portfolio volatility. Equals 1 for a single asset (or perfectly correlated
    assets) and rises as diversification lowers the portfolio vol below the
    weighted-average vol.
    """
    n = _check_cov(cov)
    sig = [math.sqrt(cov[i][i]) for i in range(n)]
    wavg = sum(weights[i] * sig[i] for i in range(n))
    pv = portfolio_variance(weights, cov)
    if pv <= 0.0:
        raise ValueError("portfolio variance must be positive")
    return wavg / math.sqrt(pv)


def max_diversification_weights(cov) -> list:
    """Most-diversified portfolio: maximizes the diversification ratio.

    The maximizer of ``(w^T sigma) / sqrt(w^T C w)`` is proportional to
    ``C^{-1} sigma`` (the tangency portfolio in the assets' own volatilities),
    normalized to sum to 1. Fully invested; may be long/short.
    """
    n = _check_cov(cov)
    inv = _invert(cov)
    sig = [math.sqrt(cov[i][i]) for i in range(n)]
    z = _matvec(inv, sig)
    total = sum(z)
    if abs(total) < 1e-14:
        raise ValueError("degenerate covariance for max diversification")
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
