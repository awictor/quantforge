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


def sobol_average_strike_rqmc(S, t, r, sigma, option_type=OptionType.CALL,
                              b=None, n_steps=6, n_paths=4096, n_rand=24,
                              seed=None) -> MCResult:
    """Randomized-QMC average-strike Asian option with an honest standard error.

    The strike is the realized arithmetic average of the monitored path, so a
    call pays ``max(S_T - A, 0)`` and a put ``max(A - S_T, 0)`` with ``A`` the
    average over the ``n_steps`` monitoring dates. Normals come from one
    ``n_steps``-dim Sobol point through the Brownian bridge, randomized by a
    per-dimension Cranley-Patterson rotation, so ``n_rand`` shifts give a genuine
    SE. The discretely-monitored analogue of
    :func:`quantforge.average_strike_asian_mc`, which it cross-checks. ``n_steps``
    is capped by the Sobol generator's dimension.
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
            avg_sum = 0.0
            s = S
            for i in range(n_steps):
                s = S * math.exp((b - 0.5 * sigma * sigma) * ((i + 1) * dt)
                                 + sigma * W[i])
                avg_sum += s
            avg = avg_sum / n_steps
            total += max(sign * (s - avg), 0.0)
        estimates.append(disc * total / n_paths)

    price, se = _summarize(estimates)
    return MCResult(price=price, std_error=se, n_paths=n_rand * n_paths)


def sobol_cliquet_rqmc(S, t, r, sigma, reset_times, local_cap=None,
                       local_floor=0.0, global_cap=None, global_floor=0.0,
                       b=None, n_paths=4096, n_rand=24, seed=None) -> MCResult:
    """Randomized-QMC capped cliquet (ratchet) with an honest standard error.

    The same product as :func:`quantforge.capped_cliquet_mc`: the payoff sums the
    per-period returns over the consecutive reset windows, each clipped to
    ``[local_floor, local_cap]``, then clips the running sum to
    ``[global_floor, global_cap]``, discounted at ``r``.

    Each path's Brownian motion at the reset times comes from one Sobol point via
    a bridge on the reset grid (:func:`_bridge_on_times`), and the per-period
    standardized shock is the bridge increment divided by ``sqrt(dt_i)``. A
    per-dimension Cranley-Patterson rotation randomizes the point set, so
    ``n_rand`` shifts give a genuine SE. The number of reset periods is capped by
    the Sobol generator's dimension. Cross-checks
    :func:`quantforge.capped_cliquet_mc`.
    """
    if b is None:
        b = r
    if sigma <= 0 or S <= 0 or t <= 0:
        raise ValueError("S, sigma, t must be positive")
    times = list(reset_times)
    n_per = len(times)
    if n_per < 1 or n_per > len(_MINIT):
        raise ValueError(f"number of reset periods must be in 1..{len(_MINIT)}")
    grid = [0.0] + times
    if any(grid[i] >= grid[i + 1] for i in range(len(grid) - 1)):
        raise ValueError("reset_times must be strictly increasing and positive")
    if abs(times[-1] - t) > 1e-9:
        raise ValueError("last reset time must be the maturity t")
    if n_rand < 2:
        raise ValueError("n_rand must be >= 2 to estimate a standard error")
    disc = math.exp(-r * t)
    rng = random.Random(seed)

    def clip(x, lo, hi):
        if lo is not None:
            x = max(x, lo)
        if hi is not None:
            x = min(x, hi)
        return x

    def payoff_from_W(W):
        total = 0.0
        w_prev = 0.0
        for i in range(n_per):
            dt = grid[i + 1] - grid[i]
            dW = W[i] - w_prev
            w_prev = W[i]
            ret = math.exp((b - 0.5 * sigma * sigma) * dt + sigma * dW) - 1.0
            total += clip(ret, local_floor, local_cap)
        return disc * clip(total, global_floor, global_cap)

    estimates = []
    for _ in range(n_rand):
        shift = [rng.random() for _ in range(n_per)]
        sob = Sobol(n_per)
        total = 0.0
        for _ in range(n_paths):
            pt = sob.next()
            u = [(pt[d] + shift[d]) % 1.0 for d in range(n_per)]
            W = _bridge_on_times(u, times)
            total += payoff_from_W(W)
        estimates.append(total / n_paths)

    price, se = _summarize(estimates)
    return MCResult(price=price, std_error=se, n_paths=n_rand * n_paths)


def sobol_double_knockout_rqmc(S, K, t, r, sigma, lower, upper,
                               option_type=OptionType.CALL, b=None, rebate=0.0,
                               n_steps=6, n_paths=4096, n_rand=24,
                               seed=None) -> MCResult:
    """Randomized-QMC double-knockout (corridor) option with an honest standard error.

    Pays the vanilla payoff only if the spot stays strictly inside
    ``(lower, upper)`` at every one of the ``n_steps`` monitoring dates; if
    either barrier is breached it knocks out and pays the cash ``rebate`` at
    expiry. Normals come from an ``n_steps``-dim Sobol point through the Brownian
    bridge, randomized by a per-dimension Cranley-Patterson rotation, so
    ``n_rand`` shifts give a genuine SE. The discretely-monitored analogue of
    :func:`quantforge.double_knockout_mc`, which it cross-checks. ``n_steps`` is
    capped by the Sobol generator's dimension.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if not (lower < S < upper):
        raise ValueError("require lower < S < upper")
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
            knocked = False
            s = S
            for i in range(n_steps):
                s = S * math.exp((b - 0.5 * sigma * sigma) * ((i + 1) * dt)
                                 + sigma * W[i])
                if s <= lower or s >= upper:
                    knocked = True
                    break
            total += rebate if knocked else max(sign * (s - K), 0.0)
        estimates.append(disc * total / n_paths)

    price, se = _summarize(estimates)
    return MCResult(price=price, std_error=se, n_paths=n_rand * n_paths)


def sobol_autocallable_rqmc(S, t, r, sigma, observation_times, autocall_barrier,
                            coupon, protection_barrier=None, notional=1.0,
                            b=None, n_paths=4096, n_rand=24,
                            seed=None) -> MCResult:
    """Randomized-QMC autocallable structured note with an honest standard error.

    The same product as :func:`quantforge.autocallable_mc`: at each observation
    date, if the spot is at or above ``autocall_barrier`` the note redeems early
    paying ``notional (1 + coupon k)`` (``k`` = observation number), discounted;
    if it never autocalls, at maturity the holder gets the notional back unless
    the spot finished below ``protection_barrier`` (a down-and-in put on the
    notional), taking ``notional S_T / S`` instead.

    Each path's ``len(observation_times)`` Brownian values come from one Sobol
    point through the Brownian bridge (dominant variance on the leading, most
    uniform coordinates), and a per-dimension Cranley-Patterson rotation
    randomizes the point set, so ``n_rand`` shifts give a genuine SE. The number
    of observations is capped by the Sobol generator's dimension. Cross-checks
    :func:`quantforge.autocallable_mc`.
    """
    if b is None:
        b = r
    obs = list(observation_times)
    n_obs = len(obs)
    if n_obs < 1 or n_obs > len(_MINIT):
        raise ValueError(f"number of observations must be in 1..{len(_MINIT)}")
    if any(obs[i] >= obs[i + 1] for i in range(n_obs - 1)) or obs[0] <= 0:
        raise ValueError("observation_times must be strictly increasing and positive")
    if abs(obs[-1] - t) > 1e-9:
        raise ValueError("last observation must be the maturity t")
    if n_rand < 2:
        raise ValueError("n_rand must be >= 2 to estimate a standard error")
    if sigma <= 0 or S <= 0 or t <= 0:
        raise ValueError("S, sigma, t must be positive")
    rng = random.Random(seed)

    # The Brownian bridge builds W at times k*dt on a uniform grid; run it on a
    # grid fine enough that every observation date lands on a node.
    def payoff_from_W(W):
        # W has length n_obs, the Brownian motion at each observation time.
        for k in range(n_obs):
            tau = obs[k]
            s = S * math.exp((b - 0.5 * sigma * sigma) * tau + sigma * W[k])
            if s >= autocall_barrier and k < n_obs - 1:
                return math.exp(-r * tau) * notional * (1.0 + coupon * (k + 1))
        disc = math.exp(-r * t)
        if protection_barrier is not None and s < protection_barrier:
            return disc * notional * (s / S)
        return disc * notional * (1.0 + coupon * n_obs)

    # Build the Brownian motion at the (possibly non-uniform) observation dates
    # from independent normals via the incremental construction, using the Sobol
    # point mapped through the bridge on the observation grid.
    estimates = []
    for _ in range(n_rand):
        shift = [rng.random() for _ in range(n_obs)]
        sob = Sobol(n_obs)
        total = 0.0
        for _ in range(n_paths):
            pt = sob.next()
            u = [(pt[d] + shift[d]) % 1.0 for d in range(n_obs)]
            # Bridge on the observation times directly.
            W = _bridge_on_times(u, obs)
            total += payoff_from_W(W)
        estimates.append(total / n_paths)

    price, se = _summarize(estimates)
    return MCResult(price=price, std_error=se, n_paths=n_rand * n_paths)


def _bridge_on_times(unifs, times):
    """Brownian motion at arbitrary increasing ``times`` from uniforms via a bridge.

    Fills the endpoint first (from the leading, most-uniform Sobol coordinate),
    then successive midpoints, so the dominant variance loads onto the early
    coordinates. Returns ``[W(times[0]), ..., W(times[-1])]`` with ``W(0)=0``.
    """
    n = len(times)
    W = [0.0] * n
    T = times[-1]
    filled = [False] * n
    order = []

    def bisect(lo, hi):
        # lo, hi are indices into times; -1 denotes t=0 (W=0).
        if hi - lo <= 1:
            return
        mid = (lo + hi) // 2
        order.append((mid, lo, hi))
        bisect(lo, mid)
        bisect(mid, hi)

    order.append((n - 1, -1, -1))       # endpoint uses coordinate 0
    bisect(-1, n - 1)

    u_idx = 0
    for (idx, lo, hi) in order:
        z = norm_ppf(unifs[u_idx])
        u_idx += 1
        if lo == -1 and hi == -1:
            W[n - 1] = math.sqrt(T) * z
        else:
            t_lo = 0.0 if lo == -1 else times[lo]
            t_hi = times[hi]
            t_mid = times[idx]
            w_lo = 0.0 if lo == -1 else W[lo]
            w_hi = W[hi]
            mean = w_lo + (w_hi - w_lo) * (t_mid - t_lo) / (t_hi - t_lo)
            var = (t_hi - t_mid) * (t_mid - t_lo) / (t_hi - t_lo)
            W[idx] = mean + math.sqrt(var) * z
    return W


def sobol_barrier_digital_rqmc(S, K, H, t, r, sigma, option_type=OptionType.CALL,
                               barrier="up-in", b=None, cash=1.0, n_steps=6,
                               n_paths=4096, n_rand=24, seed=None) -> MCResult:
    """Randomized-QMC barrier-contingent cash-or-nothing digital, honest SE.

    Pays ``cash`` at expiry iff the option finishes in the money (call
    ``S_T > K``, put ``S_T < K``) AND the barrier condition holds over the
    ``n_steps`` monitoring dates: ``up-in``/``down-in`` need the barrier touched,
    ``up-out``/``down-out`` need it untouched ("up" watches ``S >= H``, "down"
    ``S <= H``). Normals come from an ``n_steps``-dim Sobol point through the
    Brownian bridge, randomized by a per-dimension Cranley-Patterson rotation, so
    ``n_rand`` shifts give a genuine SE. The discrete analogue of
    :func:`quantforge.barrier_digital_mc`, which it cross-checks.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if H <= 0:
        raise ValueError("barrier H must be positive")
    if b is None:
        b = r
    barrier = str(barrier).lower()
    if barrier not in ("up-in", "down-in", "up-out", "down-out"):
        raise ValueError("barrier must be up-in/down-in/up-out/down-out")
    if n_steps < 1 or n_steps > len(_MINIT):
        raise ValueError(f"n_steps must be in 1..{len(_MINIT)}")
    if n_rand < 2:
        raise ValueError("n_rand must be >= 2 to estimate a standard error")
    up = barrier.startswith("up")
    knock_in = barrier.endswith("in")
    call = ot is OptionType.CALL
    dt = t / n_steps
    disc = math.exp(-r * t)
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
            touched = (up and S >= H) or (not up and S <= H)
            s = S
            for i in range(n_steps):
                s = S * math.exp((b - 0.5 * sigma * sigma) * ((i + 1) * dt)
                                 + sigma * W[i])
                if (up and s >= H) or (not up and s <= H):
                    touched = True
            barrier_ok = touched if knock_in else not touched
            itm = (s > K) if call else (s < K)
            if barrier_ok and itm:
                total += cash
        estimates.append(disc * total / n_paths)

    price, se = _summarize(estimates)
    return MCResult(price=price, std_error=se, n_paths=n_rand * n_paths)


def sobol_barrier_rqmc(S, K, H, t, r, sigma, option_type=OptionType.CALL,
                       barrier="down-out", b=None, rebate=0.0, n_steps=6,
                       n_paths=4096, n_rand=24, seed=None) -> MCResult:
    """Randomized-QMC discretely-monitored single-barrier option with honest SE.

    Prices the knock-out/knock-in vanilla barrier option monitored at the
    ``n_steps`` dates. ``barrier`` is ``down-out``/``down-in``/``up-out``/
    ``up-in``; "down" watches ``S <= H``, "up" watches ``S >= H``. ``rebate`` is
    paid at expiry to killed knock-outs or never-activated knock-ins.

    Each path's normals come from one ``n_steps``-dimensional Sobol point through
    the Brownian bridge, and a per-dimension Cranley-Patterson rotation
    randomizes the point set, so ``n_rand`` shifts give i.i.d. QMC estimates
    whose spread is a genuine SE. This is the discretely-monitored analogue of
    :func:`quantforge.barrier_mc` with ``brownian_bridge=False`` (it does not add
    the continuity correction), and cross-checks it. ``n_steps`` is capped by the
    Sobol generator's dimension.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if H <= 0:
        raise ValueError("barrier H must be positive")
    if b is None:
        b = r
    barrier = str(barrier).lower()
    if barrier not in ("down-out", "down-in", "up-out", "up-in"):
        raise ValueError("barrier must be down-out/down-in/up-out/up-in")
    if n_steps < 1 or n_steps > len(_MINIT):
        raise ValueError(f"n_steps must be in 1..{len(_MINIT)}")
    if n_rand < 2:
        raise ValueError("n_rand must be >= 2 to estimate a standard error")
    up = barrier.startswith("up")
    knock_in = barrier.endswith("in")
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
            touched = (up and S >= H) or (not up and S <= H)
            s = S
            for i in range(n_steps):
                s = S * math.exp((b - 0.5 * sigma * sigma) * ((i + 1) * dt)
                                 + sigma * W[i])
                if (up and s >= H) or (not up and s <= H):
                    touched = True
            alive = (not touched) if not knock_in else touched
            if alive:
                total += max(sign * (s - K), 0.0)
            else:
                total += rebate
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
