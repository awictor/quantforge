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
        dens.append(risk_neutral_density_from_smile(S0, t, r, vol_fn, K, q=q,
                                                    dK=max(dK, 1e-4 * S0)))
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
