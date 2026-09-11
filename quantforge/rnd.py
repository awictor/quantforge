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


def smile_arbitrage_violations(S0, t, r, vol_fn, q=0.0, n=400, width=8.0,
                               tol=1e-9):
    """Strikes where an implied-vol smile has butterfly (density) arbitrage.

    Scans a log-moneyness grid of ``width`` forward standard deviations and returns
    the strikes where the Breeden-Litzenberger density
    (:func:`risk_neutral_density_from_smile`) is negative beyond ``-tol``. A
    negative density means the call price is locally concave in strike -- a
    butterfly spread with a negative cost -- so the smile admits static arbitrage
    there. An empty list means the smile is butterfly-arbitrage-free on the grid.
    Model-free: works for any ``vol_fn`` (SVI, SABR, vanna-volga, raw quotes).
    """
    F = S0 * math.exp((r - q) * t)
    sd = vol_fn(F) * math.sqrt(t)
    lo = F * math.exp(-width * sd)
    hi = F * math.exp(width * sd)
    dK = (hi - lo) / n
    bad = []
    for i in range(n + 1):
        K = lo + i * dK
        h = min(max(1e-4 * S0, 0.25 * dK), 0.5 * K)
        g = risk_neutral_density_from_smile(S0, t, r, vol_fn, K, q=q, dK=h)
        if g < -tol:
            bad.append(K)
    return bad


def smile_is_arbitrage_free(S0, t, r, vol_fn, q=0.0, n=400, width=8.0,
                            tol=1e-9) -> bool:
    """True if the smile has no butterfly arbitrage on the scanned grid.

    Convenience wrapper: ``not smile_arbitrage_violations(...)``.
    """
    return not smile_arbitrage_violations(S0, t, r, vol_fn, q=q, n=n,
                                          width=width, tol=tol)


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


def risk_neutral_var_from_smile(S0, t, r, vol_fn, alpha=0.99, q=0.0, dK=None,
                                width=12.0):
    """Risk-neutral Value-at-Risk of the terminal simple return over ``[0, t]``.

    Works with the loss ``L = 1 - S_T / S0`` (a positive number is a loss). The
    ``alpha``-VaR is the ``alpha``-quantile of ``L``: with probability ``alpha``
    the loss does not exceed it. Since ``L <= v`` iff ``S_T >= S0 (1 - v)``,

        VaR_alpha = 1 - Q(1 - alpha) / S0,

    where ``Q`` is :func:`risk_neutral_quantile_from_smile`. Returned as a
    positive fraction of ``S0`` (e.g. ``0.18`` = an 18% loss). A flat smile
    matches the lognormal VaR.
    """
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie strictly in (0, 1)")
    K = risk_neutral_quantile_from_smile(S0, t, r, vol_fn, 1.0 - alpha, q=q,
                                         dK=dK, width=width)
    return 1.0 - K / S0


def risk_neutral_cvar_from_smile(S0, t, r, vol_fn, alpha=0.99, q=0.0, n=8000,
                                 width=12.0):
    """Risk-neutral Conditional VaR (expected shortfall) of the terminal return.

    ``CVaR_alpha = E[L | L >= VaR_alpha]``, the average loss in the worst
    ``1 - alpha`` tail. With the loss ``L = 1 - S_T / S0`` and the tail threshold
    ``K = Q(1 - alpha)`` (so ``P(S_T <= K) = 1 - alpha``),

        CVaR_alpha = 1 - E[S_T ; S_T <= K] / (S0 (1 - alpha)),

    where the truncated expectation ``E[S_T ; S_T <= K] = integral_0^K S g(S) dS``
    is taken against the Breeden-Litzenberger density ``g``. Returned as a
    positive fraction of ``S0`` and always ``>= VaR_alpha``. A flat smile matches
    the lognormal expected shortfall.
    """
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie strictly in (0, 1)")
    tail = 1.0 - alpha
    K = risk_neutral_quantile_from_smile(S0, t, r, vol_fn, tail, q=q, width=width)
    # Truncated first moment E[S_T ; S_T <= K] via trapezoid on the lower tail.
    lo = 1e-6 * S0
    dK = (K - lo) / n
    trunc = 0.0
    prev = lo * max(risk_neutral_density_from_smile(
        S0, t, r, vol_fn, lo, q=q, dK=min(1e-4 * S0, 0.5 * lo)), 0.0)
    for i in range(1, n + 1):
        S = lo + i * dK
        h = min(max(1e-4 * S0, 0.25 * dK), 0.5 * S)
        cur = S * max(risk_neutral_density_from_smile(
            S0, t, r, vol_fn, S, q=q, dK=h), 0.0)
        trunc += 0.5 * (prev + cur) * dK
        prev = cur
    return 1.0 - trunc / (S0 * tail)


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
