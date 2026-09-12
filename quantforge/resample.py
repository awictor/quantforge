"""Bootstrap and jackknife resampling for statistic confidence intervals.

Nonparametric confidence intervals for a sample statistic by resampling: the IID
bootstrap (resample observations with replacement), the stationary bootstrap of
Politis-Romano (geometric-length blocks, for serially-correlated series), and the
delete-one jackknife. All use a deterministic linear-congruential stream so
results are reproducible per seed. Pure standard library.
"""

import math


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
