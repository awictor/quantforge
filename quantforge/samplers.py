"""Random variate samplers for common distributions.

Direct samplers for the normal (Box-Muller), exponential (inverse CDF), gamma
(Marsaglia-Tsang), and Poisson (Knuth) distributions, driven by a deterministic
linear-congruential stream so runs are reproducible per seed. Each returns a list of
draws whose sample moments match the distribution's. Pure standard library.
"""

import math


class _RNG:
    def __init__(self, seed):
        self.state = seed & 0x7FFFFFFF or 1

    def uniform(self):
        self.state = (1103515245 * self.state + 12345) & 0x7FFFFFFF
        # Avoid exact 0 (log of 0) by shifting into (0, 1).
        return (self.state + 0.5) / 0x80000000


def sample_normal(n, mu=0.0, sigma=1.0, seed=1234567):
    """``n`` normal draws with mean ``mu`` and standard deviation ``sigma`` (Box-Muller)."""
    if n < 0:
        raise ValueError("n must be non-negative")
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    rng = _RNG(seed)
    out = []
    while len(out) < n:
        u1 = rng.uniform()
        u2 = rng.uniform()
        r = math.sqrt(-2.0 * math.log(u1))
        out.append(mu + sigma * r * math.cos(2.0 * math.pi * u2))
        out.append(mu + sigma * r * math.sin(2.0 * math.pi * u2))
    return out[:n]


def sample_exponential(n, rate=1.0, seed=1234567):
    """``n`` exponential draws with the given ``rate`` (mean ``1/rate``), inverse-CDF method."""
    if n < 0:
        raise ValueError("n must be non-negative")
    if rate <= 0:
        raise ValueError("rate must be positive")
    rng = _RNG(seed)
    return [-math.log(rng.uniform()) / rate for _ in range(n)]


def _one_gamma(rng, shape):
    """A single Gamma(shape, 1) draw (Marsaglia-Tsang), shape > 0."""
    if shape < 1.0:
        # Boost: Gamma(a) = Gamma(a+1) * U^(1/a).
        u = rng.uniform()
        return _one_gamma(rng, shape + 1.0) * u ** (1.0 / shape)
    d = shape - 1.0 / 3.0
    c = 1.0 / math.sqrt(9.0 * d)
    while True:
        # Standard normal via Box-Muller.
        u1 = rng.uniform()
        u2 = rng.uniform()
        x = math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)
        v = (1.0 + c * x) ** 3
        if v <= 0:
            continue
        u = rng.uniform()
        if math.log(u) < 0.5 * x * x + d - d * v + d * math.log(v):
            return d * v


def sample_gamma(n, shape, scale=1.0, seed=1234567):
    """``n`` gamma draws with the given ``shape`` (``k``) and ``scale`` (``theta``).

    Mean ``shape*scale``, variance ``shape*scale^2``. Uses Marsaglia-Tsang (with the
    small-shape boost). ``shape`` and ``scale`` must be positive.
    """
    if n < 0:
        raise ValueError("n must be non-negative")
    if shape <= 0 or scale <= 0:
        raise ValueError("shape and scale must be positive")
    rng = _RNG(seed)
    return [_one_gamma(rng, shape) * scale for _ in range(n)]


def sample_poisson(n, lam, seed=1234567):
    """``n`` Poisson draws with mean ``lam`` (Knuth's algorithm)."""
    if n < 0:
        raise ValueError("n must be non-negative")
    if lam < 0:
        raise ValueError("lam must be non-negative")
    rng = _RNG(seed)
    threshold = math.exp(-lam)
    out = []
    for _ in range(n):
        k = 0
        p = 1.0
        while True:
            p *= rng.uniform()
            if p <= threshold:
                break
            k += 1
        out.append(k)
    return out
