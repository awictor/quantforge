"""n-asset Levy moment-matched basket option.

A basket option on ``sum_i w_i S_i`` has no closed form because a sum of
lognormals is not lognormal. Levy's approximation matches the first two moments of
the basket's terminal value to a single lognormal and prices it with
Black-Scholes. Exact for one asset; a fast, accurate approximation for a
positively-weighted basket. Generalizes the two-asset ``basket_option`` to any
number of assets. Pure standard library.
"""

import math

from .mathfns import norm_cdf


def levy_basket_option(spots, weights, strike, t, r, sigmas, corr,
                       q=None, is_call=True):
    """Levy moment-matched basket call/put on ``sum_i w_i S_i``.

    Parameters
    ----------
    spots, weights, sigmas : sequences of length ``n``.
    strike, t, r : option strike, maturity, risk-free rate.
    corr : ``n x n`` correlation matrix.
    q : optional per-asset dividend yields (defaults to zeros).
    is_call : call if True, else put.

    Returns
    -------
    float
        Basket option value. Reduces to Black-Scholes for a single asset.
    """
    n = len(spots)
    if not (len(weights) == len(sigmas) == n):
        raise ValueError("spots, weights, sigmas must have equal length")
    if len(corr) != n or any(len(row) != n for row in corr):
        raise ValueError("corr must be n x n")
    if t < 0:
        raise ValueError("t must be non-negative")
    if q is None:
        q = [0.0] * n

    # Forward of each leg and the basket's first moment.
    F = [weights[i] * spots[i] * math.exp((r - q[i]) * t) for i in range(n)]
    M1 = sum(F)
    disc = math.exp(-r * t)
    if M1 <= 0.0:
        # Degenerate basket (net short); no positive-lognormal match.
        raise ValueError("basket forward must be positive")
    if t == 0.0:
        payoff = max(M1 - strike, 0.0) if is_call else max(strike - M1, 0.0)
        return disc * payoff

    # Second moment: E[B^2] = sum_ij F_i F_j exp(rho_ij sigma_i sigma_j t).
    M2 = 0.0
    for i in range(n):
        for j in range(n):
            M2 += F[i] * F[j] * math.exp(corr[i][j] * sigmas[i] * sigmas[j] * t)

    # Match to a lognormal: effective vol and forward.
    v = math.log(M2 / (M1 * M1))          # = sigma_B^2 * t
    if v <= 0.0:
        payoff = max(M1 - strike, 0.0) if is_call else max(strike - M1, 0.0)
        return disc * payoff
    sig_t = math.sqrt(v)
    d1 = (math.log(M1 / strike) + 0.5 * v) / sig_t
    d2 = d1 - sig_t
    if is_call:
        return disc * (M1 * norm_cdf(d1) - strike * norm_cdf(d2))
    return disc * (strike * norm_cdf(-d2) - M1 * norm_cdf(-d1))
