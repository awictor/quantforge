"""Risk-neutral density from an implied-vol smile, and payoff pricing.

Breeden & Litzenberger (1978): the risk-neutral density of the terminal spot is
the discounted second strike-derivative of the call price,

    g(K) = e^{r t} d^2 C / dK^2,

so any European payoff can be priced model-free from a smile by integrating it
against ``g``:

    price = e^{-r t} integral payoff(K) g(K) dK.

Here the call prices come from a smile ``vol_fn(K)`` via Black-Scholes and the
density is a central second difference on a fine strike grid. A flat smile
recovers the lognormal density and reprices vanillas/digitals exactly. Pure
standard library.
"""

import math


def risk_neutral_density_from_smile(S0, t, r, vol_fn, K, q=0.0, dK=None):
    """Breeden-Litzenberger risk-neutral density ``g(K)`` at strike ``K``.

    ``g(K) = e^{r t} d^2 C / dK^2`` with the call priced at the smile vol; a
    negative value flags butterfly arbitrage in the smile there.
    """
    from .bsm import call_price
    if dK is None:
        dK = 1e-3 * S0

    def C(k):
        return call_price(S0, k, t, r, vol_fn(k), b=r - q)

    d2 = (C(K + dK) - 2.0 * C(K) + C(K - dK)) / (dK * dK)
    return math.exp(r * t) * d2


def risk_neutral_cdf_from_smile(S0, t, r, vol_fn, K, q=0.0, dK=None):
    """Risk-neutral CDF ``F(K) = P(S_T <= K)`` implied by an implied-vol smile.

    From Breeden-Litzenberger, the digital-put price is ``e^{-rt} P(S_T <= K)`` and
    equals ``-dC/dK`` discounted, so

        F(K) = 1 + e^{r t} dC/dK,

    with the call priced at the smile vol ``vol_fn`` and ``dC/dK`` a central
    difference. Clamped to ``[0, 1]`` (a value hitting the clamp flags a smile
    that is not arbitrage-free at ``K``). A flat smile recovers the Black-Scholes
    ``N(-d2)``.
    """
    from .bsm import call_price
    if dK is None:
        dK = 1e-3 * S0
    h = min(dK, 0.5 * K)

    def C(k):
        return call_price(S0, k, t, r, vol_fn(k), b=r - q)

    dCdK = (C(K + h) - C(K - h)) / (2.0 * h)
    F = 1.0 + math.exp(r * t) * dCdK
    return min(1.0, max(0.0, F))


def risk_neutral_quantile_from_smile(S0, t, r, vol_fn, p, q=0.0, dK=None,
                                     width=12.0, tol=1e-8, max_iter=200):
    """Inverse risk-neutral CDF: the strike ``K`` with ``P(S_T <= K) = p``.

    Bisection on :func:`risk_neutral_cdf_from_smile` over a log-moneyness bracket
    of ``+/- width`` forward standard deviations. ``p`` in ``(0, 1)``. Requires the
    smile CDF to be monotone on the bracket (true for an arbitrage-free smile).
    """
    if not 0.0 < p < 1.0:
        raise ValueError("p must lie strictly in (0, 1)")
    F = S0 * math.exp((r - q) * t)
    sd = vol_fn(F) * math.sqrt(t)
    lo = F * math.exp(-width * sd)
    hi = F * math.exp(width * sd)

    def cdf(K):
        return risk_neutral_cdf_from_smile(S0, t, r, vol_fn, K, q=q, dK=dK)

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        fm = cdf(mid)
        if abs(fm - p) < tol or (hi - lo) < tol * F:
            return mid
        if fm < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def density_grid_from_smile(S0, t, r, vol_fn, q=0.0, n=400, width=8.0):
    """Return ``(strikes, density)`` of the risk-neutral density on a grid.

    Strikes span ``width`` standard deviations of log-moneyness around the
    forward. Useful for plotting or integrating custom payoffs.
    """
    F = S0 * math.exp((r - q) * t)
    sd = vol_fn(F) * math.sqrt(t)
    lo = F * math.exp(-width * sd)
    hi = F * math.exp(width * sd)
    dK = (hi - lo) / n
    strikes, dens = [], []
    for i in range(n + 1):
        K = lo + i * dK
        strikes.append(K)
        # The FD step must stay strictly inside (0, K) so K - h > 0 at the low
        # nodes (a wide grid can have the grid spacing exceed the first strike).
        h = min(max(1e-4 * S0, 0.25 * dK), 0.5 * K)
        dens.append(risk_neutral_density_from_smile(S0, t, r, vol_fn, K, q=q,
                                                    dK=h))
    return strikes, dens


def price_payoff_from_density(S0, t, r, vol_fn, payoff, q=0.0, n=400, width=8.0):
    """Price a European payoff ``payoff(S_T)`` model-free from the smile density.

    ``price = e^{-r t} integral payoff(K) g(K) dK`` by the trapezoidal rule over
    the Breeden-Litzenberger density grid. A flat smile reprices vanillas and
    digitals to Black-Scholes.
    """
    strikes, dens = density_grid_from_smile(S0, t, r, vol_fn, q=q, n=n,
                                            width=width)
    disc = math.exp(-r * t)
    total = 0.0
    for i in range(len(strikes) - 1):
        k0, k1 = strikes[i], strikes[i + 1]
        g0 = payoff(k0) * max(dens[i], 0.0)
        g1 = payoff(k1) * max(dens[i + 1], 0.0)
        total += 0.5 * (g0 + g1) * (k1 - k0)
    return disc * total
