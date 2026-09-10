"""Sobol low-discrepancy sequence and a Brownian-bridge QMC path pricer.

Sobol points fill the unit cube far more evenly than pseudo-random draws,
especially in the leading dimensions, giving close to ``O((log N)^d / N)`` error.
For a *path*-dependent payoff the dimension is the number of time steps, and the
crucial trick is the **Brownian bridge** construction: it loads most of each
path's variance onto the first few Sobol coordinates (the endpoint, then the
midpoint, then quarter points, ...), where Sobol is most uniform, instead of
spreading it evenly across a random-walk increment per coordinate.

This module implements a Sobol generator from primitive-polynomial direction
numbers (Joe-Kuo initial values for the first dimensions) via the Gray-code
recurrence, and a Brownian-bridge Asian/European Monte Carlo that draws its
normals from Sobol points. Pure standard library.
"""

import math
import random
from typing import List

from .mathfns import norm_ppf
from .bsm import OptionType, _coerce_type, _validate
from .montecarlo import MCResult, _summarize

_BITS = 30
_SCALE = float(1 << _BITS)

# Primitive polynomials (as the integer of their inner coefficient bits, a_i)
# and Joe-Kuo initial direction numbers m_i for dimensions 2.. (dim 1 is trivial).
# Dimension 1 uses v_k = 2^{BITS-1-k}; dims below cover the leading coordinates
# that carry the Brownian bridge's dominant variance.
_POLY = [0, 1, 1, 2, 1, 4]                        # a-bits per dim (dim>=2)
_MINIT = [
    [],                    # dim 1 (special-cased)
    [1],                   # dim 2, poly x+1        (degree 1)
    [1, 1],                # dim 3, poly x^2+x+1    (degree 2)
    [1, 3, 7],             # dim 4, poly x^3+x+1
    [1, 1, 5],             # dim 5, poly x^3+x^2+1
    [1, 3, 1, 1],          # dim 6, poly x^4+x+1
]


def _direction_numbers(dim: int) -> List[List[int]]:
    """Direction-number table v[d][k] (scaled to 2^BITS) for ``dim`` dimensions."""
    V = []
    for d in range(dim):
        v = [0] * _BITS
        if d == 0:
            for k in range(_BITS):
                v[k] = 1 << (_BITS - 1 - k)
        else:
            m = _MINIT[d][:]           # initial m_1..m_s
            s = len(m)
            a = _POLY[d]
            for k in range(s):
                v[k] = m[k] << (_BITS - 1 - k)
            for k in range(s, _BITS):
                val = v[k - s] ^ (v[k - s] >> s)
                for j in range(1, s):
                    if (a >> (s - 1 - j)) & 1:
                        val ^= v[k - j]
                v[k] = val
        V.append(v)
    return V


class Sobol:
    """A Sobol sequence generator (Gray-code recurrence)."""

    def __init__(self, dim: int):
        if dim < 1 or dim > len(_MINIT):
            raise ValueError(f"dim must be in 1..{len(_MINIT)}")
        self.dim = dim
        self._V = _direction_numbers(dim)
        self._x = [0] * dim
        self._count = 0

    def next(self) -> List[float]:
        """Return the next Sobol point in [0, 1)^dim.

        The first call returns the origin's successor (index 1), skipping the
        degenerate all-zeros point, which is the usual convention.
        """
        c = self._count
        self._count += 1
        if c == 0:
            return [0.5] * self.dim   # index 1 point; avoids the origin
        # Index of the lowest zero bit of c (Gray-code recurrence).
        j = 0
        cc = c
        while cc & 1:
            cc >>= 1
            j += 1
        out = []
        for d in range(self.dim):
            self._x[d] ^= self._V[d][j]
            out.append(self._x[d] / _SCALE)
        return out


def _bridge_order(n: int):
    """Brownian-bridge fill order and (left, right, this) index triples.

    Returns a list of ``(idx, left, right)`` where ``left``/``right`` are already
    -set neighbouring time indices (or -1 for the known zero at t=0), applied in
    an order that resolves the endpoint first, then successive midpoints.
    """
    order = [(n - 1, -1, -1)]           # endpoint first (uses dimension 0)
    # Recursively bisect [ -1 (t0), n-1 ].
    stack = [(-1, n - 1)]
    while stack:
        lo, hi = stack.pop(0)
        if hi - lo <= 1:
            continue
        mid = (lo + hi) // 2
        order.append((mid, lo, hi))
        stack.append((lo, mid))
        stack.append((mid, hi))
    return order


def brownian_bridge_path(unifs: List[float], t: float):
    """Build a Brownian path W_0..W_n at times k*dt from Sobol uniforms.

    ``unifs`` has length ``n`` (one per time step). The endpoint is drawn from
    the first coordinate, then successive midpoints from the rest -- so the
    leading (most uniform) Sobol dimensions carry the dominant variance.
    Returns the list ``[W_1, ..., W_n]`` (W_0 = 0).
    """
    n = len(unifs)
    dt = t / n
    W = [0.0] * (n + 1)   # W[0] = 0 at t=0
    order = _bridge_order(n)
    u_idx = 0
    for (idx, lo, hi) in order:
        z = norm_ppf(unifs[u_idx])
        u_idx += 1
        if lo == -1 and hi == -1:
            # Endpoint W[n] ~ N(0, t).
            W[n] = math.sqrt(t) * z
        else:
            t_lo = 0.0 if lo == -1 else (lo + 1) * dt
            t_hi = (hi + 1) * dt
            t_mid = (idx + 1) * dt
            w_lo = 0.0 if lo == -1 else W[lo + 1]
            w_hi = W[hi + 1]
            # Bridge: mean interpolates, variance is the conditional bridge var.
            mean = w_lo + (w_hi - w_lo) * (t_mid - t_lo) / (t_hi - t_lo)
            var = (t_hi - t_mid) * (t_mid - t_lo) / (t_hi - t_lo)
            W[idx + 1] = mean + math.sqrt(var) * z
    return W[1:]


def sobol_european(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                   n_paths=8192):
    """QMC European price using a 1-D Sobol sequence (single-step payoff).

    A direct Sobol analogue of :func:`quantforge.european_qmc`; converges faster
    than pseudo-random Monte Carlo for this smooth one-dimensional integral.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    drift = (b - 0.5 * sigma * sigma) * t
    vol = sigma * math.sqrt(t)
    disc = math.exp(-r * t)
    sign = 1.0 if ot is OptionType.CALL else -1.0
    sob = Sobol(1)
    total = 0.0
    for _ in range(n_paths):
        u = sob.next()[0]
        sT = S * math.exp(drift + vol * norm_ppf(u))
        total += max(sign * (sT - K), 0.0)
    return disc * total / n_paths


def sobol_european_rqmc(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                        n_paths=4096, n_rand=24, seed=None) -> MCResult:
    """Randomized-QMC European price with an honest standard error.

    Plain Sobol QMC (:func:`sobol_european`) returns a single number with no
    error estimate -- the points are deterministic, so there is no variance to
    report. Randomized QMC restores an unbiased error bar by applying a
    Cranley-Patterson rotation: shift the whole Sobol point set by a random
    ``U ~ Uniform[0,1)`` modulo 1. Each shift preserves the sequence's low
    discrepancy but makes the resulting estimate an unbiased draw, so ``n_rand``
    independent shifts give ``n_rand`` i.i.d. QMC estimates whose spread is a
    genuine standard error.

    The returned :class:`MCResult` has ``price`` = mean over the randomizations,
    ``std_error`` = their across-randomization SE, and ``n_paths`` = the total
    points evaluated (``n_rand * n_paths``). For this smooth 1-D integral the
    RQMC SE falls off far faster than pseudo-random Monte Carlo's ``1/sqrt(N)``.

    Cross-checks the closed-form Black-Scholes value.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if n_rand < 2:
        raise ValueError("n_rand must be >= 2 to estimate a standard error")
    drift = (b - 0.5 * sigma * sigma) * t
    vol = sigma * math.sqrt(t)
    disc = math.exp(-r * t)
    sign = 1.0 if ot is OptionType.CALL else -1.0
    rng = random.Random(seed)

    estimates = []
    for _ in range(n_rand):
        shift = rng.random()
        sob = Sobol(1)
        total = 0.0
        for _ in range(n_paths):
            u = (sob.next()[0] + shift) % 1.0        # Cranley-Patterson rotation
            sT = S * math.exp(drift + vol * norm_ppf(u))
            total += max(sign * (sT - K), 0.0)
        estimates.append(disc * total / n_paths)

    price, se = _summarize(estimates)
    return MCResult(price=price, std_error=se, n_paths=n_rand * n_paths)


def sobol_asian_rqmc(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                     n_steps=6, n_paths=4096, n_rand=24, seed=None) -> MCResult:
    """Randomized-QMC arithmetic Asian price with an honest standard error.

    The multi-dimensional analogue of :func:`sobol_european_rqmc`. Each path's
    ``n_steps`` normals come from one ``n_steps``-dimensional Sobol point mapped
    through the Brownian bridge (dominant variance on the leading, most uniform
    coordinates). A per-dimension Cranley-Patterson rotation -- shift each Sobol
    coordinate by an independent ``U ~ Uniform[0,1)`` modulo 1 -- randomizes the
    whole point set without disturbing its low discrepancy, so ``n_rand``
    independent shifts give i.i.d. QMC estimates whose spread is a genuine SE.

    Returns an :class:`MCResult` with the mean price, the across-randomization
    SE, and ``n_paths`` = total points (``n_rand * n_paths``). Cross-checks the
    geometric-control-variate :func:`quantforge.arithmetic_asian_mc`.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if n_steps < 1 or n_steps > len(_MINIT):
        raise ValueError(f"n_steps must be in 1..{len(_MINIT)}")
    if n_rand < 2:
        raise ValueError("n_rand must be >= 2 to estimate a standard error")
    dt = t / n_steps
    disc = math.exp(-r * t)
    sign = 1.0 if ot is OptionType.CALL else -1.0
    rng = random.Random(seed)

    estimates = []
    for _ in range(n_rand):
        shift = [rng.random() for _ in range(n_steps)]
        sob = Sobol(n_steps)
        total = 0.0
        for _ in range(n_paths):
            pt = sob.next()
            u = [(pt[d] + shift[d]) % 1.0 for d in range(n_steps)]
            W = brownian_bridge_path(u, t)
            avg = 0.0
            for i in range(n_steps):
                tk = (i + 1) * dt
                s = S * math.exp((b - 0.5 * sigma * sigma) * tk + sigma * W[i])
                avg += s
            avg /= n_steps
            total += max(sign * (avg - K), 0.0)
        estimates.append(disc * total / n_paths)

    price, se = _summarize(estimates)
    return MCResult(price=price, std_error=se, n_paths=n_rand * n_paths)


def sobol_lookback_rqmc(S, t, r, sigma, option_type=OptionType.CALL, b=None,
                        n_steps=6, n_paths=4096, n_rand=24,
                        seed=None) -> MCResult:
    """Randomized-QMC floating-strike lookback with an honest standard error.

    Prices the discretely-monitored floating-strike lookback -- call payoff
    ``S_T - min_i S_{t_i}``, put payoff ``max_i S_{t_i} - S_T`` -- where the
    running extreme is taken over the ``n_steps`` monitoring dates (plus the
    known ``S_0``). Each path's normals come from one ``n_steps``-dimensional
    Sobol point through the Brownian bridge, and a per-dimension
    Cranley-Patterson rotation randomizes the point set, so ``n_rand`` shifts
    give i.i.d. QMC estimates whose spread is a genuine SE.

    Discrete monitoring always *under*-prices the continuously-monitored
    Goldman-Sosin-Gatto :func:`quantforge.floating_strike_lookback` (fewer
    sampling dates see less extreme highs/lows); the gap shrinks as ``n_steps``
    grows. ``n_steps`` is capped by the Sobol generator's dimension. Returns an
    :class:`MCResult` with the mean price, across-randomization SE, and
    ``n_paths`` = total points.
    """
    ot = _coerce_type(option_type)
    _validate(S, S, t, sigma)
    if b is None:
        b = r
    if n_steps < 1 or n_steps > len(_MINIT):
        raise ValueError(f"n_steps must be in 1..{len(_MINIT)}")
    if n_rand < 2:
        raise ValueError("n_rand must be >= 2 to estimate a standard error")
    dt = t / n_steps
    disc = math.exp(-r * t)
    call = ot is OptionType.CALL
    rng = random.Random(seed)

    estimates = []
    for _ in range(n_rand):
        shift = [rng.random() for _ in range(n_steps)]
        sob = Sobol(n_steps)
        total = 0.0
        for _ in range(n_paths):
            pt = sob.next()
            u = [(pt[d] + shift[d]) % 1.0 for d in range(n_steps)]
            W = brownian_bridge_path(u, t)
            smin = smax = S
            sT = S
            for i in range(n_steps):
                tk = (i + 1) * dt
                sT = S * math.exp((b - 0.5 * sigma * sigma) * tk + sigma * W[i])
                if sT < smin:
                    smin = sT
                if sT > smax:
                    smax = sT
            total += (sT - smin) if call else (smax - sT)
        estimates.append(disc * total / n_paths)

    price, se = _summarize(estimates)
    return MCResult(price=price, std_error=se, n_paths=n_rand * n_paths)


def sobol_fixed_lookback_rqmc(S, K, t, r, sigma, option_type=OptionType.CALL,
                              b=None, n_steps=6, n_paths=4096, n_rand=24,
                              seed=None) -> MCResult:
    """Randomized-QMC discrete fixed-strike lookback with an honest standard error.

    Prices the discretely-monitored fixed-strike lookback -- call payoff
    ``max(max_i S_{t_i} - K, 0)``, put payoff ``max(K - min_i S_{t_i}, 0)`` --
    where the running extreme is taken over ``S_0`` and the ``n_steps``
    monitoring dates. Normals come from an ``n_steps``-dimensional Sobol point
    through the Brownian bridge, randomized by a per-dimension Cranley-Patterson
    rotation so ``n_rand`` shifts give a genuine SE.

    Discrete monitoring *under*-prices the continuously-monitored Conze-
    Viswanathan :func:`quantforge.fixed_strike_lookback` (fewer sampling dates
    see less extreme highs/lows); the gap shrinks as ``n_steps`` grows.
    ``n_steps`` is capped by the Sobol generator's dimension.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if n_steps < 1 or n_steps > len(_MINIT):
        raise ValueError(f"n_steps must be in 1..{len(_MINIT)}")
    if n_rand < 2:
        raise ValueError("n_rand must be >= 2 to estimate a standard error")
    dt = t / n_steps
    disc = math.exp(-r * t)
    call = ot is OptionType.CALL
    rng = random.Random(seed)

    estimates = []
    for _ in range(n_rand):
        shift = [rng.random() for _ in range(n_steps)]
        sob = Sobol(n_steps)
        total = 0.0
        for _ in range(n_paths):
            pt = sob.next()
            u = [(pt[d] + shift[d]) % 1.0 for d in range(n_steps)]
            W = brownian_bridge_path(u, t)
            smin = smax = S
            for i in range(n_steps):
                tk = (i + 1) * dt
                s = S * math.exp((b - 0.5 * sigma * sigma) * tk + sigma * W[i])
                if s < smin:
                    smin = s
                if s > smax:
                    smax = s
            total += max(smax - K, 0.0) if call else max(K - smin, 0.0)
        estimates.append(disc * total / n_paths)

    price, se = _summarize(estimates)
    return MCResult(price=price, std_error=se, n_paths=n_rand * n_paths)


def sobol_asian(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                n_steps=6, n_paths=8192):
    """QMC arithmetic-average Asian price with Sobol + a Brownian bridge.

    Each path's ``n_steps`` normals come from one Sobol point mapped through the
    Brownian-bridge construction, so the dominant path variance lands on the
    leading (most uniform) Sobol dimensions. ``n_steps`` is capped by the
    generator's dimension.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if n_steps < 1 or n_steps > len(_MINIT):
        raise ValueError(f"n_steps must be in 1..{len(_MINIT)}")
    dt = t / n_steps
    disc = math.exp(-r * t)
    sign = 1.0 if ot is OptionType.CALL else -1.0
    sob = Sobol(n_steps)
    total = 0.0
    for _ in range(n_paths):
        u = sob.next()
        W = brownian_bridge_path(u, t)   # Brownian motion at each step
        avg = 0.0
        for i in range(n_steps):
            tk = (i + 1) * dt
            s = S * math.exp((b - 0.5 * sigma * sigma) * tk + sigma * W[i])
            avg += s
        avg /= n_steps
        total += max(sign * (avg - K), 0.0)
    return disc * total / n_paths
