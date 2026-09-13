"""Sharpe's returns-based style analysis (RBSA).

Explains a fund's returns as a portfolio of asset-class indices, with weights that
are non-negative and sum to one (a long-only, fully-invested mix). It solves the
constrained least-squares problem

    min_w  Var(fund - sum_k w_k index_k)   s.t.  w_k >= 0,  sum_k w_k = 1,

by projected gradient descent on the simplex. The resulting weights are the fund's
implied style exposures, and the unexplained variance is its selection (active)
return. Pure standard library.
"""


def _project_simplex(v):
    """Euclidean projection of ``v`` onto the probability simplex (Duchi et al.)."""
    n = len(v)
    u = sorted(v, reverse=True)
    css = 0.0
    rho = 0
    theta = 0.0
    for i in range(n):
        css += u[i]
        t = (css - 1.0) / (i + 1)
        if u[i] - t > 0:
            rho = i + 1
            theta = t
    return [max(x - theta, 0.0) for x in v]


def style_analysis(fund_returns, index_returns, max_iter=5000, lr=None):
    """Returns-based style analysis: implied long-only index weights of a fund.

    ``fund_returns`` is the return series; ``index_returns`` is a list of index
    return series (one per style factor), each aligned with the fund. Returns a dict
    with ``weights`` (non-negative, summing to one), ``r_squared`` (fraction of fund
    variance explained by the style mix), and ``tracking_error`` (stdev of the
    unexplained residual). Solved by projected-gradient descent on the simplex.
    """
    n = len(fund_returns)
    k = len(index_returns)
    if k == 0:
        raise ValueError("need at least one index")
    if any(len(idx) != n for idx in index_returns):
        raise ValueError("all series must be the same length as the fund")
    if n < 2:
        raise ValueError("need at least two observations")

    # Gradient of the residual sum of squares in w: 2 X'(X w - y).
    # X is n x k (index returns), y is the fund.
    def residual_ss(w):
        ss = 0.0
        for t in range(n):
            pred = sum(w[j] * index_returns[j][t] for j in range(k))
            ss += (fund_returns[t] - pred) ** 2
        return ss

    # Step size from the largest index second moment (Lipschitz-ish).
    if lr is None:
        max_norm = max(sum(index_returns[j][t] ** 2 for t in range(n))
                       for j in range(k))
        lr = 1.0 / (2.0 * max_norm) if max_norm > 0 else 0.1

    w = [1.0 / k] * k
    for _ in range(max_iter):
        grad = [0.0] * k
        for t in range(n):
            resid = sum(w[j] * index_returns[j][t] for j in range(k)) - fund_returns[t]
            for j in range(k):
                grad[j] += 2.0 * resid * index_returns[j][t]
        w_new = _project_simplex([w[j] - lr * grad[j] for j in range(k)])
        if max(abs(w_new[j] - w[j]) for j in range(k)) < 1e-12:
            w = w_new
            break
        w = w_new

    # Diagnostics.
    resid = [fund_returns[t] - sum(w[j] * index_returns[j][t] for j in range(k))
             for t in range(n)]
    mean_f = sum(fund_returns) / n
    ss_tot = sum((fund_returns[t] - mean_f) ** 2 for t in range(n))
    ss_res = sum(r * r for r in resid)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    mean_r = sum(resid) / n
    te = (sum((r - mean_r) ** 2 for r in resid) / (n - 1)) ** 0.5
    return {"weights": w, "r_squared": r_squared, "tracking_error": te}
