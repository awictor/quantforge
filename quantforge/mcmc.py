"""Markov chain Monte Carlo samplers: random-walk Metropolis and Hamiltonian Monte Carlo.

MCMC draws samples from a target density known only up to a constant -- the workhorse of
Bayesian inference. Two samplers are provided:

* :func:`metropolis_hastings` -- random-walk Metropolis. Propose a Gaussian perturbation, accept
  with probability ``min(1, exp(logp_new - logp_old))``. Simple and robust; mixes slowly in high
  dimensions.
* :func:`hamiltonian_monte_carlo` -- HMC. Introduces a momentum variable and simulates
  Hamiltonian dynamics with leapfrog steps, using the *gradient* of the log-density to make long,
  low-rejection moves. The gradient comes from reverse-mode autodiff, so the log-density is
  written once with :class:`quantforge.reverse_ad.Var` arithmetic.

Both use a seeded :class:`quantforge.pcg.PCG32` stream for reproducibility. Pure standard
library.
"""

import math

from .pcg import PCG32
from .particle_filter import pcg_gaussian
from .reverse_jacobian import reverse_gradient_vector


def metropolis_hastings(log_prob, x0, n_samples, step=0.5, seed=12345, burn_in=0, thin=1):
    """Random-walk Metropolis sampler for a target ``log_prob``.

    ``log_prob`` maps a length-``d`` list of floats to the log target density (up to an additive
    constant). Proposals are isotropic Gaussian with standard deviation ``step``. Returns a dict
    with ``samples`` (list of accepted states after burn-in/thinning) and ``accept_rate``.
    """
    rng = PCG32(seed)
    d = len(x0)
    x = [float(v) for v in x0]
    lp = log_prob(x)
    samples = []
    n_accept = 0
    total = burn_in + n_samples * thin
    for it in range(total):
        prop = [x[i] + pcg_gaussian(rng, 0.0, step) for i in range(d)]
        lp_prop = log_prob(prop)
        if math.log(rng.random() + 1e-300) < lp_prop - lp:
            x, lp = prop, lp_prop
            n_accept += 1
        if it >= burn_in and (it - burn_in) % thin == 0:
            samples.append(x[:])
    return {"samples": samples, "accept_rate": n_accept / total}


def hamiltonian_monte_carlo(log_prob, x0, n_samples, step=0.1, n_leapfrog=20,
                            seed=12345, burn_in=0):
    """Hamiltonian Monte Carlo sampler using autodiff gradients of ``log_prob``.

    ``log_prob`` maps a length-``d`` list of :class:`quantforge.reverse_ad.Var` to a single
    ``Var`` (the log target density up to a constant). Each iteration samples a Gaussian momentum,
    runs ``n_leapfrog`` leapfrog steps of step size ``step`` on the Hamiltonian
    ``H = -log_prob(x) + |p|^2 / 2``, and Metropolis-accepts on the total energy. Returns a dict
    with ``samples`` and ``accept_rate``.
    """
    rng = PCG32(seed)
    d = len(x0)

    def value(z):
        from .reverse_ad import Var
        return log_prob([Var(zi) for zi in z]).value

    def grad(z):
        return reverse_gradient_vector(log_prob, z)   # d(log_prob)/dz

    x = [float(v) for v in x0]
    samples = []
    n_accept = 0
    total = burn_in + n_samples
    for it in range(total):
        p = [pcg_gaussian(rng, 0.0, 1.0) for _ in range(d)]
        x_new = x[:]
        g = grad(x_new)
        p_new = [p[i] + 0.5 * step * g[i] for i in range(d)]   # half kick (grad of +logp)
        for _ in range(n_leapfrog):
            x_new = [x_new[i] + step * p_new[i] for i in range(d)]   # drift
            g = grad(x_new)
            p_new = [p_new[i] + step * g[i] for i in range(d)]       # full kick
        # correct the last half kick
        p_new = [p_new[i] - 0.5 * step * g[i] for i in range(d)]

        # Metropolis on total energy (potential U = -logp, kinetic K = |p|^2/2)
        cur_U = -value(x)
        new_U = -value(x_new)
        cur_K = 0.5 * sum(pi * pi for pi in p)
        new_K = 0.5 * sum(pi * pi for pi in p_new)
        if math.log(rng.random() + 1e-300) < (cur_U + cur_K) - (new_U + new_K):
            x = x_new
            n_accept += 1
        if it >= burn_in:
            samples.append(x[:])
    return {"samples": samples, "accept_rate": n_accept / total}


def _mean(samples, j):
    return sum(s[j] for s in samples) / len(samples)


def sample_mean(samples):
    """Component-wise mean of a list of sample vectors."""
    d = len(samples[0])
    return [_mean(samples, j) for j in range(d)]


def sample_cov(samples):
    """Component-wise covariance matrix of a list of sample vectors (population form)."""
    n = len(samples)
    d = len(samples[0])
    mu = sample_mean(samples)
    cov = [[0.0] * d for _ in range(d)]
    for s in samples:
        for i in range(d):
            for j in range(d):
                cov[i][j] += (s[i] - mu[i]) * (s[j] - mu[j])
    return [[cov[i][j] / n for j in range(d)] for i in range(d)]
