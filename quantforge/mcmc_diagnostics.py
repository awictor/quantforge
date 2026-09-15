"""Convergence diagnostics for MCMC output.

A finished MCMC run (:mod:`quantforge.mcmc`) is only trustworthy if the chain has *converged* to
the target and the samples are not too correlated. Three standard diagnostics:

* :func:`gelman_rubin` -- the R-hat statistic. Run several independent chains; if they have
  converged to the same distribution, the between-chain variance matches the within-chain
  variance and ``R-hat -> 1``. Values above ~1.1 signal non-convergence.
* :func:`integrated_autocorrelation_time` -- ``tau = 1 + 2 sum_k rho_k``, the number of steps
  between effectively independent samples, using Geyer's initial-positive-sequence truncation.
* :func:`effective_sample_size` -- ``n / tau``, how many independent draws the correlated chain
  is worth.

Pure standard library (reuses the autocovariance idea behind :mod:`quantforge.acf`).
"""


def _mean(x):
    return sum(x) / len(x)


def _autocov(x, k, mean):
    n = len(x)
    return sum((x[t] - mean) * (x[t + k] - mean) for t in range(n - k)) / n


def gelman_rubin(chains):
    """Gelman-Rubin potential-scale-reduction ``R-hat`` for a list of scalar chains.

    ``chains`` is a list of ``m >= 2`` equal-length sequences (one per independent run). Returns
    ``R-hat = sqrt(V_hat / W)`` where ``W`` is the mean within-chain variance and ``V_hat`` the
    variance estimate mixing in the between-chain spread. ``-> 1`` at convergence.
    """
    m = len(chains)
    if m < 2:
        raise ValueError("need at least 2 chains")
    n = len(chains[0])
    if n < 2 or any(len(c) != n for c in chains):
        raise ValueError("chains must have equal length >= 2")
    means = [_mean(c) for c in chains]
    grand = _mean(means)
    # within-chain variance W (unbiased per chain, averaged)
    def var(c, mu):
        return sum((v - mu) ** 2 for v in c) / (n - 1)
    W = sum(var(chains[j], means[j]) for j in range(m)) / m
    # between-chain variance B
    B = n * sum((means[j] - grand) ** 2 for j in range(m)) / (m - 1)
    if W == 0.0:
        return 1.0
    var_hat = (n - 1) / n * W + B / n
    return (var_hat / W) ** 0.5


def integrated_autocorrelation_time(x, max_lag=None):
    """Integrated autocorrelation time ``tau = 1 + 2 sum_k rho_k`` (Geyer initial positive seq).

    Sums the autocorrelations, truncating at the first lag where consecutive-pair sums turn
    negative (Geyer's rule), which keeps the estimate stable. ``tau >= 1``; larger means more
    correlated samples. Returns ``1.0`` for a zero-variance series.
    """
    n = len(x)
    if n < 2:
        raise ValueError("need at least 2 points")
    mean = _mean(x)
    g0 = _autocov(x, 0, mean)
    if g0 == 0.0:
        return 1.0
    if max_lag is None:
        # cap the work: autocorrelation of a usable chain decays well before this
        max_lag = min(n - 1, max(1000, int(10 * n ** 0.5)))
    # Geyer's initial positive sequence: sum pairs (rho_{2k} + rho_{2k+1}) while positive.
    # Autocovariances are computed lazily so we stop as soon as the pair sum turns negative
    # (typically after a handful of lags), instead of materializing all n-1 lags.
    def rho(k):
        return _autocov(x, k, mean) / g0
    tau = 1.0
    k = 1
    while k + 1 <= max_lag:
        pair = rho(k) + rho(k + 1)
        if pair <= 0:
            break
        tau += 2.0 * pair
        k += 2
    return tau


def effective_sample_size(x, max_lag=None):
    """Effective sample size ``n / tau`` -- independent-draw equivalent of a correlated chain.

    Uses :func:`integrated_autocorrelation_time`. Equals ``n`` for white noise and shrinks as
    the chain's autocorrelation grows. Never exceeds ``n``.
    """
    n = len(x)
    tau = integrated_autocorrelation_time(x, max_lag)
    return min(n, n / tau)
