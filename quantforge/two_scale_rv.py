"""Two-scale realized variance (Zhang-Mykland-Ait-Sahalia 2005).

Observed high-frequency prices are the efficient (log) price plus i.i.d.
microstructure noise: ``p_i = p*_i + eps_i``. The naive realized variance
``sum (p_{i+1} - p_i)^2`` is then badly *upward*-biased -- each squared return
carries ``2 * Var(eps)`` of pure noise, and that bias grows with the sampling
frequency, so faster sampling makes it worse, not better.

The two-scale estimator removes the bias by combining two time scales:

  * a *slow* scale that averages the realized variance over ``K`` subsampled
    grids (each taking every ``K``-th observation), which is nearly noise-free but
    uses fewer points;
  * a *fast* scale (the all-observations RV) that is noisy but whose noise bias is
    estimable, ``2 * n * Var(eps)`` with ``Var(eps)`` recovered as ``RV_fast /
    (2 n)``.

    TSRV = RV_avg - (nbar / n) * RV_fast,

with a small-sample factor ``1 / (1 - nbar/n)`` applied. This is consistent for
the integrated variance even under noise. Pure standard library.
"""


def _rv(prices):
    return sum((prices[i + 1] - prices[i]) ** 2 for i in range(len(prices) - 1))


def realized_variance_naive(prices):
    """Naive all-observations realized variance of a (log) price series.

    Upward-biased by ``2 n Var(eps)`` under i.i.d. microstructure noise.
    """
    if len(prices) < 2:
        raise ValueError("need at least 2 prices")
    return _rv(prices)


def noise_variance_estimate(prices):
    """Estimate the microstructure-noise variance ``Var(eps)``.

    Under the noise model the fast-scale RV is dominated by noise, so
    ``Var(eps) ~= RV_fast / (2 * (n - 1))`` where ``n - 1`` is the number of
    returns. Consistent as the noise dwarfs the signal at the finest scale.
    """
    if len(prices) < 2:
        raise ValueError("need at least 2 prices")
    m = len(prices) - 1
    return _rv(prices) / (2.0 * m)


def two_scale_realized_variance(prices, K=None):
    """Two-scale realized variance (TSRV), robust to microstructure noise.

    Parameters
    ----------
    prices : sequence of float
        Observed log-prices on a fine grid.
    K : int, optional
        Number of subsampling grids for the slow scale. Defaults to
        ``max(2, round(n ** (1/3)))``, the rate-optimal choice.

    Returns
    -------
    float
        Bias-corrected estimate of the integrated variance. On noise-free data it
        essentially reproduces the realized variance; under noise it is far less
        biased than :func:`realized_variance_naive`.
    """
    n = len(prices)
    if n < 3:
        raise ValueError("need at least 3 prices")
    if K is None:
        K = max(2, round(n ** (1.0 / 3.0)))
    if not (1 < K < n):
        raise ValueError("K must satisfy 1 < K < n")

    # Slow scale: average RV over the K subsampled grids (stride K, offset g).
    slow_sum = 0.0
    grids = 0
    for g in range(K):
        sub = prices[g::K]
        if len(sub) >= 2:
            slow_sum += _rv(sub)
            grids += 1
    rv_avg = slow_sum / grids

    # Fast scale: all-observations RV, and its per-grid return count.
    rv_fast = _rv(prices)
    n_returns = n - 1
    nbar = (n - K + 1) / K          # average #returns per subsampled grid

    tsrv = rv_avg - (nbar / n_returns) * rv_fast
    # Small-sample adjustment (Zhang-Mykland-Ait-Sahalia).
    adj = 1.0 / (1.0 - nbar / n_returns)
    return adj * tsrv
