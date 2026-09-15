"""Adaptive MCMC: Haario adaptive Metropolis and the slice sampler.

Plain random-walk Metropolis (:func:`quantforge.mcmc.metropolis_hastings`) needs a hand-tuned
proposal scale, and a single isotropic scale mixes badly when the target is strongly correlated
or has very different marginal scales. Two samplers here remove that tuning burden:

* :func:`adaptive_metropolis` -- Haario, Saksman & Tamminen (2001). The proposal covariance is
  continuously re-estimated from the chain's own history (``2.38^2 / d`` times the empirical
  covariance, the asymptotically optimal scaling), so it learns the target's shape on the fly.
* :func:`slice_sample` -- Neal's (2003) univariate slice sampler with stepping-out and shrinkage.
  It has no proposal scale to tune at all: it samples uniformly under the density curve.

Both use a seeded :class:`quantforge.pcg.PCG32` stream. Pure standard library.
"""

import math

from .pcg import PCG32
from .particle_filter import pcg_gaussian
from .linalg import cholesky


def adaptive_metropolis(log_prob, x0, n_samples, seed=12345, burn_in=1000,
                        init_scale=0.1, adapt_start=200, epsilon=1e-6):
    """Haario adaptive-Metropolis sampler for a target ``log_prob`` (list of floats -> logp).

    The proposal is Gaussian with covariance ``(2.38^2 / d) * Cov(history) + epsilon I`` once at
    least ``adapt_start`` samples have accumulated; before that it is isotropic with standard
    deviation ``init_scale``. Returns a dict with ``samples`` (after ``burn_in``) and
    ``accept_rate``.
    """
    rng = PCG32(seed)
    d = len(x0)
    x = [float(v) for v in x0]
    lp = log_prob(x)
    sd = (2.38 ** 2) / d
    chol = None                      # cached Cholesky factor of the current proposal cov
    samples = []
    n_accept = 0
    total = burn_in + n_samples

    # running mean and unnormalized covariance (M2) updated incrementally, O(d^2)/step,
    # so adaptation does not cost O(history) each iteration.
    count = 1
    run_mean = x[:]
    M2 = [[0.0] * d for _ in range(d)]

    for it in range(total):
        if it >= adapt_start and count > 1:
            C = [[M2[i][j] / (count - 1) for j in range(d)] for i in range(d)]
            for i in range(d):
                C[i][i] += epsilon
            prop_cov = [[sd * C[i][j] for j in range(d)] for i in range(d)]
            try:
                chol = cholesky(prop_cov)
            except (ValueError, ZeroDivisionError):
                chol = None
        if chol is not None:
            z = [pcg_gaussian(rng, 0.0, 1.0) for _ in range(d)]
            step = [sum(chol[i][k] * z[k] for k in range(i + 1)) for i in range(d)]
            prop = [x[i] + step[i] for i in range(d)]
        else:
            prop = [x[i] + pcg_gaussian(rng, 0.0, init_scale) for i in range(d)]
        lp_prop = log_prob(prop)
        if math.log(rng.random() + 1e-300) < lp_prop - lp:
            x, lp = prop, lp_prop
            n_accept += 1
        # incremental Welford update of running mean / M2 with the new state x
        count += 1
        delta = [x[i] - run_mean[i] for i in range(d)]
        for i in range(d):
            run_mean[i] += delta[i] / count
        delta2 = [x[i] - run_mean[i] for i in range(d)]
        for i in range(d):
            for j in range(d):
                M2[i][j] += delta[i] * delta2[j]
        if it >= burn_in:
            samples.append(x[:])
    return {"samples": samples, "accept_rate": n_accept / total}


def slice_sample(log_prob, x0, n_samples, w=1.0, seed=12345, burn_in=0, max_steps=50):
    """Univariate slice sampler (Neal 2003) for a scalar target ``log_prob(x) -> logp``.

    Draws an auxiliary height under the density, steps out an interval of initial width ``w``
    (up to ``max_steps`` expansions each side), then samples uniformly from the interval,
    shrinking on rejection. No proposal scale to tune. Returns a list of ``n_samples`` scalars.
    """
    rng = PCG32(seed)
    x = float(x0)
    out = []
    total = burn_in + n_samples
    for it in range(total):
        # auxiliary level: log(y) = log_prob(x) - Exponential(1)
        log_y = log_prob(x) - (-math.log(rng.random() + 1e-300))
        # step out
        u = rng.random()
        left = x - w * u
        right = left + w
        j = int(max_steps * rng.random())
        k = max_steps - 1 - j
        while j > 0 and log_prob(left) > log_y:
            left -= w
            j -= 1
        while k > 0 and log_prob(right) > log_y:
            right += w
            k -= 1
        # shrink
        for _ in range(100):
            x_new = left + rng.random() * (right - left)
            if log_prob(x_new) > log_y:
                x = x_new
                break
            if x_new < x:
                left = x_new
            else:
                right = x_new
        if it >= burn_in:
            out.append(x)
    return out
