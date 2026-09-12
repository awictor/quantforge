"""Bootstrap and jackknife resampling for statistic confidence intervals.

Nonparametric confidence intervals for a sample statistic by resampling: the IID
bootstrap (resample observations with replacement), the stationary bootstrap of
Politis-Romano (geometric-length blocks, for serially-correlated series), and the
delete-one jackknife. All use a deterministic linear-congruential stream so
results are reproducible per seed. Pure standard library.
"""

import math

from .mathfns import norm_cdf, norm_ppf


def _lcg(seed):
    state = seed & 0xFFFFFFFF

    def _next():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return state / 0x80000000

    return _next


def _percentile(sorted_vals, p):
    if not sorted_vals:
        raise ValueError("empty sample")
    if p <= 0:
        return sorted_vals[0]
    if p >= 100:
        return sorted_vals[-1]
    idx = (p / 100.0) * (len(sorted_vals) - 1)
    lo = int(math.floor(idx))
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = idx - lo
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[hi] * frac


def _mean(xs):
    return sum(xs) / len(xs)


def bootstrap_ci(data, statistic=None, n_boot=2000, confidence=0.95,
                 seed=1234567):
    """IID bootstrap confidence interval for a sample ``statistic``.

    Resamples ``data`` with replacement ``n_boot`` times, applies ``statistic`` to
    each resample, and returns ``(lower, point, upper)`` -- the percentile-method
    interval at ``confidence`` plus the statistic on the original sample. The
    interval brackets the point estimate and narrows as the sample grows.
    ``statistic`` defaults to the sample mean.
    """
    if statistic is None:
        statistic = _mean
    n = len(data)
    if n == 0:
        raise ValueError("data must be non-empty")
    if not (0.0 < confidence < 1.0):
        raise ValueError("confidence must be in (0, 1)")
    rand = _lcg(seed)
    reps = []
    for _ in range(n_boot):
        sample = [data[int(rand() * n)] for _ in range(n)]
        reps.append(statistic(sample))
    reps.sort()
    alpha = (1.0 - confidence) / 2.0
    lo = _percentile(reps, 100.0 * alpha)
    hi = _percentile(reps, 100.0 * (1.0 - alpha))
    return lo, statistic(data), hi


def stationary_bootstrap_ci(data, statistic=None, mean_block=10, n_boot=2000,
                            confidence=0.95, seed=1234567):
    """Stationary (Politis-Romano) bootstrap CI for serially-correlated data.

    Resamples geometric-length blocks (expected length ``mean_block``) wrapping
    around the series, preserving short-range dependence, then takes the
    percentile interval. Wider than the IID :func:`bootstrap_ci` for positively
    autocorrelated series (it does not spuriously shrink the variance).
    ``statistic`` defaults to the sample mean.
    """
    if statistic is None:
        statistic = _mean
    n = len(data)
    if n == 0:
        raise ValueError("data must be non-empty")
    if mean_block < 1:
        raise ValueError("mean_block must be >= 1")
    if not (0.0 < confidence < 1.0):
        raise ValueError("confidence must be in (0, 1)")
    p_new = 1.0 / mean_block
    rand = _lcg(seed)
    reps = []
    for _ in range(n_boot):
        sample = []
        idx = int(rand() * n)
        while len(sample) < n:
            sample.append(data[idx])
            if rand() < p_new:
                idx = int(rand() * n)   # start a new block
            else:
                idx = (idx + 1) % n     # continue the current block
        reps.append(statistic(sample))
    reps.sort()
    alpha = (1.0 - confidence) / 2.0
    lo = _percentile(reps, 100.0 * alpha)
    hi = _percentile(reps, 100.0 * (1.0 - alpha))
    return lo, statistic(data), hi


def moving_block_bootstrap_ci(data, statistic=None, block=10, n_boot=2000,
                              confidence=0.95, seed=1234567):
    """Moving-block (Kunsch) bootstrap CI for serially-correlated data.

    Resamples fixed-length overlapping blocks of length ``block`` from the series
    (wrapping at the end) and concatenates ceil(n / block) of them, truncated to
    ``n``, preserving within-block dependence. Like the stationary bootstrap it
    gives valid intervals for autocorrelated data -- wider than the IID
    :func:`bootstrap_ci` for a positively autocorrelated mean, and it agrees with
    the IID interval when ``block = 1``. ``statistic`` defaults to the sample mean.
    """
    if statistic is None:
        statistic = _mean
    n = len(data)
    if n == 0:
        raise ValueError("data must be non-empty")
    if block < 1 or block > n:
        raise ValueError("block must satisfy 1 <= block <= len(data)")
    if not (0.0 < confidence < 1.0):
        raise ValueError("confidence must be in (0, 1)")
    n_blocks = (n + block - 1) // block
    rand = _lcg(seed)
    reps = []
    for _ in range(n_boot):
        sample = []
        for _ in range(n_blocks):
            start = int(rand() * n)
            for k in range(block):
                sample.append(data[(start + k) % n])
        reps.append(statistic(sample[:n]))
    reps.sort()
    alpha = (1.0 - confidence) / 2.0
    lo = _percentile(reps, 100.0 * alpha)
    hi = _percentile(reps, 100.0 * (1.0 - alpha))
    return lo, statistic(data), hi


def bca_bootstrap_ci(data, statistic=None, n_boot=2000, confidence=0.95,
                     seed=1234567):
    """Bias-corrected accelerated (BCa) bootstrap confidence interval.

    Efron's BCa improves on the percentile method by correcting for median bias
    (``z0``, from the fraction of bootstrap replicates below the point estimate)
    and skewness (``a``, the acceleration from the jackknife). The percentiles are
    shifted:

        alpha1 = Phi(z0 + (z0 + z_lo)/(1 - a(z0 + z_lo)))
        alpha2 = Phi(z0 + (z0 + z_hi)/(1 - a(z0 + z_hi))).

    Reduces to the plain :func:`bootstrap_ci` when ``z0`` and ``a`` are zero
    (symmetric, unbiased statistic). Returns ``(lower, point, upper)``.
    """
    if statistic is None:
        statistic = _mean
    n = len(data)
    if n < 2:
        raise ValueError("need at least two observations")
    if not (0.0 < confidence < 1.0):
        raise ValueError("confidence must be in (0, 1)")
    rand = _lcg(seed)
    theta_hat = statistic(data)
    reps = []
    for _ in range(n_boot):
        sample = [data[int(rand() * n)] for _ in range(n)]
        reps.append(statistic(sample))
    reps_sorted = sorted(reps)
    # Bias correction z0 from the fraction of replicates below the estimate.
    n_below = sum(1 for r in reps if r < theta_hat)
    prop = n_below / n_boot
    prop = min(max(prop, 1.0 / (2 * n_boot)), 1.0 - 1.0 / (2 * n_boot))
    z0 = norm_ppf(prop)
    # Acceleration a from the jackknife influence values.
    thetas = []
    for i in range(n):
        loo = data[:i] + data[i + 1:]
        thetas.append(statistic(loo))
    theta_bar = sum(thetas) / n
    num = sum((theta_bar - t) ** 3 for t in thetas)
    den = 6.0 * (sum((theta_bar - t) ** 2 for t in thetas)) ** 1.5
    a = num / den if den != 0.0 else 0.0
    alpha = (1.0 - confidence) / 2.0
    z_lo = norm_ppf(alpha)
    z_hi = norm_ppf(1.0 - alpha)

    def _adjust(z):
        denom = 1.0 - a * (z0 + z)
        return norm_cdf(z0 + (z0 + z) / denom) if denom != 0.0 else norm_cdf(z0 + z)

    lo = _percentile(reps_sorted, 100.0 * _adjust(z_lo))
    hi = _percentile(reps_sorted, 100.0 * _adjust(z_hi))
    return lo, theta_hat, hi


def jackknife_estimate(data, statistic=None):
    """Delete-one jackknife estimate and standard error of a ``statistic``.

    Recomputes the statistic on each leave-one-out subsample. Returns
    ``(estimate, standard_error)`` with the bias-aware jackknife SE
    ``sqrt((n-1)/n * sum (theta_i - theta_bar)^2)``. ``statistic`` defaults to the
    sample mean.
    """
    if statistic is None:
        statistic = _mean
    n = len(data)
    if n < 2:
        raise ValueError("need at least two observations")
    thetas = []
    for i in range(n):
        loo = data[:i] + data[i + 1:]
        thetas.append(statistic(loo))
    theta_bar = sum(thetas) / n
    var = (n - 1) / n * sum((t - theta_bar) ** 2 for t in thetas)
    return statistic(data), math.sqrt(var)
