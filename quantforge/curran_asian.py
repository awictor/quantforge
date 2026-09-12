"""Curran (1994) arithmetic-average Asian approximation.

Curran's method conditions the arithmetic average on the geometric mean ``G``
(which is exactly lognormal) and integrates analytically, giving a sharper price
than the two-moment Turnbull-Wakeman match -- especially for many monitoring
dates. For ``n`` equally spaced driftless (forward-measure) monitoring dates
``t_i = i T / n`` with each leg's forward equal to ``forward``:

    ln S_i ~ N(mu_i, v_i),  mu_i = ln F - 0.5 sigma^2 t_i,  v_i = sigma^2 t_i,
    mu = mean_i mu_i,  v_x = Var(ln G) = (sigma^2 / n^2) sum_ij min(t_i, t_j),
    sigma_{iG} = Cov(ln S_i, ln G) = (sigma^2 / n) sum_j min(t_i, t_j).

A threshold ``K*`` on ``G`` is solved from the average-equals-strike condition and
the call is a sum of ``N(.)`` terms. Reduces to Black at ``n = 1``. Pure standard
library.
"""

import math

from .mathfns import norm_cdf


def curran_asian(forward, strike, sigma, r, expiry, n_avg, is_call=True):
    """Curran arithmetic-average Asian option price.

    Parameters mirror :func:`~quantforge.turnbull_wakeman_asian`. Returns the
    discounted option value; call and put satisfy
    ``C - P = e^{-rT}(forward - strike)``.
    """
    if forward <= 0 or strike <= 0:
        raise ValueError("forward and strike must be positive")
    if expiry < 0 or sigma < 0:
        raise ValueError("expiry and sigma must be non-negative")
    if n_avg < 1:
        raise ValueError("n_avg must be a positive integer")
    disc = math.exp(-r * expiry)
    n = n_avg

    if expiry == 0.0 or sigma == 0.0:
        payoff = max(forward - strike, 0.0) if is_call else max(strike - forward, 0.0)
        return disc * payoff

    dt = expiry / n
    t = [(i + 1) * dt for i in range(n)]
    mu_i = [math.log(forward) - 0.5 * sigma * sigma * t[i] for i in range(n)]
    v_i = [sigma * sigma * t[i] for i in range(n)]
    mu = sum(mu_i) / n
    v_x = 0.0
    for i in range(n):
        for j in range(n):
            v_x += min(t[i], t[j])
    v_x *= sigma * sigma / (n * n)
    sig_iG = [sigma * sigma / n * sum(min(t[i], t[j]) for j in range(n))
              for i in range(n)]
    sx = math.sqrt(v_x)

    # Solve K*: (1/n) sum_i exp(mu_i + (sig_iG/v_x)(lnK* - mu) + 0.5(v_i - sig_iG^2/v_x)) = K
    def avg_given_lnG(lnK):
        s = 0.0
        for i in range(n):
            s += math.exp(mu_i[i] + (sig_iG[i] / v_x) * (lnK - mu)
                          + 0.5 * (v_i[i] - sig_iG[i] * sig_iG[i] / v_x))
        return s / n

    lo, hi = math.log(strike) - 10.0 * sx - 5.0, math.log(strike) + 10.0 * sx + 5.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if avg_given_lnG(mid) < strike:
            lo = mid
        else:
            hi = mid
    lnKstar = 0.5 * (lo + hi)

    d = (mu - lnKstar) / sx
    call = 0.0
    for i in range(n):
        di = d + sig_iG[i] / sx
        call += math.exp(mu_i[i] + 0.5 * v_i[i]) * norm_cdf(di)
    call = disc * (call / n - strike * norm_cdf(d))
    if is_call:
        return max(call, 0.0)
    return call - disc * (forward - strike)
