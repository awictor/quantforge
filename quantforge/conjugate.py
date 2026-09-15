"""Conjugate Bayesian updating for the standard one-parameter models.

When the prior and likelihood are conjugate, the posterior is available in closed form -- no
sampling needed. These are the three textbook cases plus a generic highest-density-interval
routine:

* :func:`beta_binomial_posterior` -- Beta prior + Binomial data -> Beta posterior (a proportion).
* :func:`gamma_poisson_posterior` -- Gamma prior + Poisson counts -> Gamma posterior (a rate).
* :func:`normal_normal_posterior` -- Normal prior + Normal data (known variance) -> Normal
  posterior (a mean).
* :func:`hpd_interval` -- the shortest interval containing a given posterior mass, found on a
  grid of a density callable (works for any of the above, or an arbitrary 1-D posterior).

Each posterior is returned as a dict of its (hyper)parameters plus ``mean`` and ``var``. Pure
standard library.
"""

import math


def beta_binomial_posterior(prior_alpha, prior_beta, successes, trials):
    """Beta-Binomial conjugate update for a success probability ``p``.

    Prior ``p ~ Beta(prior_alpha, prior_beta)`` and ``successes`` out of ``trials`` Bernoulli
    outcomes give posterior ``Beta(alpha + s, beta + (n - s))``. Returns the posterior ``alpha``,
    ``beta``, ``mean`` and ``var``.
    """
    if prior_alpha <= 0 or prior_beta <= 0:
        raise ValueError("prior parameters must be positive")
    if not (0 <= successes <= trials):
        raise ValueError("need 0 <= successes <= trials")
    a = prior_alpha + successes
    b = prior_beta + (trials - successes)
    mean = a / (a + b)
    var = a * b / ((a + b) ** 2 * (a + b + 1))
    return {"alpha": a, "beta": b, "mean": mean, "var": var}


def gamma_poisson_posterior(prior_shape, prior_rate, total_count, n_obs):
    """Gamma-Poisson conjugate update for a Poisson rate ``lambda``.

    Prior ``lambda ~ Gamma(shape, rate)`` (rate = inverse scale) and ``n_obs`` observations with
    summed count ``total_count`` give posterior ``Gamma(shape + total_count, rate + n_obs)``.
    Returns posterior ``shape``, ``rate``, ``mean`` and ``var``.
    """
    if prior_shape <= 0 or prior_rate <= 0:
        raise ValueError("prior parameters must be positive")
    if total_count < 0 or n_obs < 0:
        raise ValueError("counts must be non-negative")
    shape = prior_shape + total_count
    rate = prior_rate + n_obs
    mean = shape / rate
    var = shape / (rate * rate)
    return {"shape": shape, "rate": rate, "mean": mean, "var": var}


def normal_normal_posterior(prior_mean, prior_var, data, data_var):
    """Normal-Normal conjugate update for a mean ``mu`` with known observation variance.

    Prior ``mu ~ Normal(prior_mean, prior_var)`` and ``data`` (list) drawn with known variance
    ``data_var`` give a Normal posterior. Precisions add: ``1/post_var = 1/prior_var + n/data_var``
    and the posterior mean is the precision-weighted average of prior mean and sample mean.
    Returns posterior ``mean`` and ``var`` (``var`` is the variance of ``mu``, not of the data).
    """
    if prior_var <= 0 or data_var <= 0:
        raise ValueError("variances must be positive")
    n = len(data)
    if n == 0:
        return {"mean": prior_mean, "var": prior_var}
    sample_mean = sum(data) / n
    prior_prec = 1.0 / prior_var
    data_prec = n / data_var
    post_prec = prior_prec + data_prec
    post_var = 1.0 / post_prec
    post_mean = (prior_prec * prior_mean + data_prec * sample_mean) / post_prec
    return {"mean": post_mean, "var": post_var}


def hpd_interval(density, lo, hi, mass=0.95, n_grid=10000):
    """Highest-posterior-density interval of a 1-D ``density`` over ``[lo, hi]``.

    Evaluates ``density`` on a uniform grid, then finds the shortest interval (by lowering a
    horizontal threshold on the density until the enclosed probability reaches ``mass``). Returns
    ``(low, high)``. The density need not be normalized. Suits any of the conjugate posteriors
    (pass the corresponding pdf) or an arbitrary unimodal 1-D posterior.
    """
    if not 0 < mass < 1:
        raise ValueError("mass must be in (0, 1)")
    dx = (hi - lo) / n_grid
    xs = [lo + (i + 0.5) * dx for i in range(n_grid)]
    dens = [max(density(x), 0.0) for x in xs]
    total = sum(dens) * dx
    if total <= 0:
        raise ValueError("density integrates to zero on [lo, hi]")
    # Lower a threshold; keep grid points with density >= threshold until mass is covered.
    order = sorted(range(n_grid), key=lambda i: dens[i], reverse=True)
    acc = 0.0
    chosen = []
    target = mass * total
    for i in order:
        acc += dens[i] * dx
        chosen.append(xs[i])
        if acc >= target:
            break
    return (min(chosen), max(chosen))
