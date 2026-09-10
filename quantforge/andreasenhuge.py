"""Andreasen-Huge (2011) single-step arbitrage-free local-volatility smile.

Andreasen & Huge showed that one *implicit* finite-difference step of Dupire's
forward equation, from expiry back to today, produces a full set of call prices
that is arbitrage-free by construction -- monotone decreasing and convex in
strike -- for any positive local-volatility function. In the forward measure
(work with forward calls; discount outside) the undiscounted call prices
``c_i = C(K_i)`` solve the tridiagonal system

    (1 + 2 a_i) c_i - a_i c_{i-1} - a_i c_{i+1} = payoff_i,
    a_i = 0.5 * T * sigma_i^2 * K_i^2 / dx^2,   payoff_i = max(F - K_i, 0),

with Dirichlet ends ``c = F - K`` at the low-strike edge and ``c = 0`` at the
high-strike edge. The system matrix is a strictly diagonally-dominant M-matrix,
so the solution is automatically a valid (no-arbitrage) call surface. Inverting
the map -- choosing the piecewise ``sigma_i`` that reprices a market smile -- is
a fast per-strike bootstrap. Pure standard library.
"""

import math

from .bsm import OptionType, _coerce_type, call_price
from .implied import implied_volatility


def _solve_tridiag(lower, diag, upper, rhs):
    """Thomas algorithm for a tridiagonal system; returns the solution list."""
    n = len(diag)
    c = [0.0] * n
    d = [0.0] * n
    c[0] = upper[0] / diag[0]
    d[0] = rhs[0] / diag[0]
    for i in range(1, n):
        m = diag[i] - lower[i] * c[i - 1]
        c[i] = upper[i] / m if i < n - 1 else 0.0
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / m
    x = [0.0] * n
    x[n - 1] = d[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = d[i] - c[i] * x[i + 1]
    return x


def andreasen_huge_prices(F, strikes, T, local_vols):
    """Arbitrage-free forward call prices from one implicit Dupire step.

    ``strikes`` is an increasing grid, ``local_vols`` the per-strike local vol
    (same length). Returns the list of undiscounted forward call prices
    ``C(K_i)`` at expiry ``T``. Multiply by ``e^{-rT}`` for the discounted price
    if the forward already embeds the carry.

    The prices are monotone decreasing and convex in strike by construction, for
    any positive ``local_vols`` -- the whole point of the scheme.
    """
    n = len(strikes)
    if n < 3 or len(local_vols) != n:
        raise ValueError("need >=3 strikes and matching local_vols")
    payoff = [max(F - K, 0.0) for K in strikes]

    lower = [0.0] * n
    diag = [1.0] * n
    upper = [0.0] * n
    rhs = list(payoff)
    # Dirichlet ends: c[0] = F - K[0] (deep ITM), c[-1] = 0 (deep OTM).
    diag[0] = 1.0
    rhs[0] = max(F - strikes[0], 0.0)
    diag[n - 1] = 1.0
    rhs[n - 1] = 0.0
    for i in range(1, n - 1):
        dx = 0.5 * (strikes[i + 1] - strikes[i - 1])
        a = 0.5 * T * local_vols[i] ** 2 * strikes[i] ** 2 / (dx * dx)
        lower[i] = -a
        diag[i] = 1.0 + 2.0 * a
        upper[i] = -a
    return _solve_tridiag(lower, diag, upper, rhs)


def andreasen_huge_smile(F, strikes, T, local_vols, r=0.0):
    """Implied-vol smile from the Andreasen-Huge arbitrage-free call prices.

    Prices the grid with :func:`andreasen_huge_prices` (spot ``S = F e^{-rT}``,
    carry ``b = r`` so the pricing forward is ``F``) and inverts each to a
    Black-Scholes implied vol. Returns ``(log_moneyness, vol)`` pairs sorted by
    strike on the forward ``F``.
    """
    prices = andreasen_huge_prices(F, strikes, T, local_vols)
    S = F * math.exp(-r * T)
    disc = math.exp(-r * T)
    out = []
    for K, c_fwd in zip(strikes, prices):
        c = disc * c_fwd
        try:
            iv = implied_volatility(c, S, K, T, r, OptionType.CALL, b=r)
        except ValueError:
            continue
        if iv <= 1e-8:
            continue   # deep ITM/OTM edge (Dirichlet) inverts degenerately
        out.append((math.log(K / F), iv))
    return out


def andreasen_huge_calibrate(F, strikes, T, market_vols, r=0.0,
                             max_iter=60, tol=1e-8):
    """Calibrate per-strike local vols so the AH smile matches market vols.

    Bootstraps the local vol at each interior strike by a 1-D bisection on the
    single-strike implied-vol error (the implicit step couples neighbours only
    weakly, so a few sweeps converge). Returns ``(local_vols, rmse)`` with
    ``rmse`` the root-mean-square implied-vol error over the interior strikes.
    """
    n = len(strikes)
    if n < 3 or len(market_vols) != n:
        raise ValueError("need >=3 strikes and matching market_vols")
    # Seed local vols at the market implied vols.
    lv = list(market_vols)

    for _ in range(max_iter):
        smile = andreasen_huge_smile(F, strikes, T, lv, r=r)
        model = {round(k, 12): iv for k, iv in
                 [(math.log(strikes[i] / F), v) for i, v in
                  enumerate([iv for _, iv in smile])]} if False else None
        # Re-price and compute per-strike model vols aligned to strikes.
        prices = andreasen_huge_prices(F, strikes, T, lv)
        S = F * math.exp(-r * T)
        disc = math.exp(-r * T)
        worst = 0.0
        for i in range(1, n - 1):
            try:
                iv = implied_volatility(disc * prices[i], S, strikes[i], T,
                                        r, OptionType.CALL, b=r)
            except ValueError:
                continue
            err = iv - market_vols[i]
            worst = max(worst, abs(err))
            # Local vol and implied vol move together; nudge proportionally.
            lv[i] *= market_vols[i] / iv if iv > 1e-8 else 1.0
        if worst < tol:
            break

    # Final RMSE over interior strikes.
    prices = andreasen_huge_prices(F, strikes, T, lv)
    S = F * math.exp(-r * T)
    disc = math.exp(-r * T)
    sse = 0.0
    cnt = 0
    for i in range(1, n - 1):
        try:
            iv = implied_volatility(disc * prices[i], S, strikes[i], T, r,
                                    OptionType.CALL, b=r)
        except ValueError:
            continue
        sse += (iv - market_vols[i]) ** 2
        cnt += 1
    return lv, math.sqrt(sse / cnt) if cnt else float("nan")
