"""n-asset rainbow options (best-of / worst-of) by Monte Carlo.

Options on the maximum or minimum of ``n`` correlated assets at expiry:

    best-of call  = E[max(0, max_i S_i(T) - K)],
    worst-of call = E[max(0, min_i S_i(T) - K)],

both discounted. Terminal prices are simulated as correlated geometric Brownian
motions (Cholesky of the correlation matrix on standard normals); a deterministic
linear-congruential stream keeps the estimate reproducible per seed. Always
``worst-of <= single-asset <= best-of``. Pure standard library.
"""

import math

from .linalg import cholesky


def _lcg_normals(seed):
    state = seed & 0x7FFFFFFF
    cache = []

    def _uni():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return (state + 0.5) / 0x80000000

    def _next():
        if cache:
            return cache.pop()
        u1 = _uni()
        u2 = _uni()
        rr = math.sqrt(-2.0 * math.log(u1))
        cache.append(rr * math.sin(2.0 * math.pi * u2))
        return rr * math.cos(2.0 * math.pi * u2)

    return _next


def rainbow_option_mc(spots, strike, t, r, sigmas, corr, q=None,
                      best=True, is_call=True, n_paths=100000, seed=1234567):
    """Monte Carlo price of an ``n``-asset best-of / worst-of option.

    Parameters
    ----------
    spots, sigmas : per-asset spot and volatility (length ``n``).
    strike, t, r : option strike, maturity, risk-free rate.
    corr : ``n x n`` correlation matrix.
    q : optional per-asset dividend yields.
    best : True for best-of (max), False for worst-of (min).
    is_call : call if True, else put.
    n_paths, seed : simulation controls.

    Returns
    -------
    float
        Discounted Monte Carlo option value.
    """
    n = len(spots)
    if not (len(sigmas) == n):
        raise ValueError("spots and sigmas must have equal length")
    if len(corr) != n or any(len(row) != n for row in corr):
        raise ValueError("corr must be n x n")
    if t < 0:
        raise ValueError("t must be non-negative")
    if q is None:
        q = [0.0] * n
    L = cholesky(corr)                       # raises if not positive-definite
    gen = _lcg_normals(seed)
    disc = math.exp(-r * t)
    sqrt_t = math.sqrt(t)
    drift = [(r - q[i] - 0.5 * sigmas[i] * sigmas[i]) * t for i in range(n)]

    acc = 0.0
    for _ in range(n_paths):
        z = [gen() for _ in range(n)]
        term = []
        for i in range(n):
            cz = sum(L[i][j] * z[j] for j in range(i + 1))
            term.append(spots[i] * math.exp(drift[i] + sigmas[i] * sqrt_t * cz))
        chosen = max(term) if best else min(term)
        if is_call:
            acc += max(chosen - strike, 0.0)
        else:
            acc += max(strike - chosen, 0.0)
    return disc * acc / n_paths
