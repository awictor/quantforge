"""Bootstrap particle filter (sequential importance resampling).

The Kalman family (linear / extended / unscented) assumes Gaussian noise and unimodal
posteriors. A particle filter makes neither assumption: it represents the state distribution by
a weighted cloud of samples ("particles"), propagates each through the (possibly nonlinear,
non-Gaussian) dynamics, reweights them by the observation likelihood, and resamples to focus
computation where the probability mass is. This is the tool for multimodal tracking, bearings-
only problems, and heavy-tailed noise.

This is the classic bootstrap filter (Gordon-Salmond-Smith 1993): the proposal is the state
transition itself, so weights update by the measurement likelihood alone. Systematic resampling
(low-variance) is triggered when the effective sample size drops below half. A seeded PCG32
stream (:class:`quantforge.pcg.PCG32`) makes runs reproducible. Pure standard library.
"""

import math

from .pcg import PCG32


def particle_filter(observations, transition, log_likelihood, init_sampler,
                    n_particles=1000, seed=12345, resample_threshold=0.5):
    """Bootstrap particle filter over ``observations``.

    Callbacks (all receive the PCG32 ``rng`` where randomness is needed):

    * ``init_sampler(rng)`` -> an initial state (any object the other callbacks understand).
    * ``transition(state, rng)`` -> the next state, sampled from the process model.
    * ``log_likelihood(observation, state)`` -> log ``p(obs | state)``.

    Returns a dict with ``means`` (the weighted-mean state estimate at each step; states must
    support scalar or per-component averaging -- floats or equal-length lists), ``ess`` (the
    effective sample size at each step), and ``n_resample`` (how many steps resampled).
    """
    rng = PCG32(seed)
    particles = [init_sampler(rng) for _ in range(n_particles)]
    w = [1.0 / n_particles] * n_particles

    means = []
    ess_hist = []
    n_resample = 0

    for z in observations:
        # propagate through the transition model
        particles = [transition(p, rng) for p in particles]
        # reweight by the observation likelihood (log-domain, then normalize stably)
        logw = [math.log(w[i]) + log_likelihood(z, particles[i]) for i in range(n_particles)]
        mx = max(logw)
        w = [math.exp(lw - mx) for lw in logw]
        total = sum(w)
        w = [wi / total for wi in w]

        means.append(_weighted_mean(particles, w))
        ess = 1.0 / sum(wi * wi for wi in w)
        ess_hist.append(ess)

        if ess < resample_threshold * n_particles:
            particles = _systematic_resample(particles, w, rng)
            w = [1.0 / n_particles] * n_particles
            n_resample += 1

    return {"means": means, "ess": ess_hist, "n_resample": n_resample}


def _weighted_mean(particles, w):
    first = particles[0]
    if isinstance(first, (list, tuple)):
        d = len(first)
        return [sum(w[i] * particles[i][j] for i in range(len(particles))) for j in range(d)]
    return sum(w[i] * particles[i] for i in range(len(particles)))


def _systematic_resample(particles, w, rng):
    n = len(particles)
    positions = [(i + rng.random()) / n for i in range(n)]
    cumulative = []
    c = 0.0
    for wi in w:
        c += wi
        cumulative.append(c)
    cumulative[-1] = 1.0  # guard against round-off
    out = []
    j = 0
    for pos in positions:
        while pos > cumulative[j]:
            j += 1
        out.append(particles[j])
    return out


def pcg_gaussian(rng, mu=0.0, sigma=1.0):
    """Draw a normal sample from a PCG32 stream via the Box-Muller transform.

    Convenience for writing ``transition``/``init_sampler`` callbacks without pulling in another
    RNG. Uses two uniforms; returns a single normal deviate.
    """
    u1 = rng.random()
    u2 = rng.random()
    if u1 < 1e-300:
        u1 = 1e-300
    return mu + sigma * math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)
