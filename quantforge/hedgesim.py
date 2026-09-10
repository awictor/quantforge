"""Discrete delta-hedging P&L simulator.

Sells one European option and hedges it by holding ``delta`` shares of the
underlying, rebalancing at ``n_steps`` equally spaced dates along a simulated
GBM path. In the Black-Scholes limit of continuous rehedging the terminal P&L
is exactly zero (the option premium replicates the payoff); with discrete
rehedging a residual **hedging error** remains, and its standard deviation
scales like ``1/sqrt(n_steps)`` (Boyle-Emanuel). This module quantifies that
error distribution by Monte Carlo.

Convention: we are **short** the option (collect the premium) and **long**
delta shares. Cash earns the risk-free rate between rebalances. Terminal P&L =
final cash + share value - option payoff owed.
"""

import math
import random
from dataclasses import dataclass

from .bsm import price as bsm_price, delta as bsm_delta, OptionType, _coerce_type, _validate


@dataclass(frozen=True)
class HedgeResult:
    mean_pnl: float
    std_pnl: float          # the hedging-error standard deviation
    min_pnl: float
    max_pnl: float
    n_paths: int
    n_steps: int

    def percentile(self, samples, p):  # pragma: no cover - convenience only
        s = sorted(samples)
        i = min(len(s) - 1, max(0, int(p * len(s))))
        return s[i]


def simulate_delta_hedge(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                         n_steps=50, n_paths=20_000, hedge_vol=None,
                         real_vol=None, seed=None, return_samples=False):
    """Monte Carlo a discretely delta-hedged short option position.

    Args:
        n_steps: rebalancing dates over the option's life.
        hedge_vol: volatility used to compute the hedge delta (defaults to
            ``sigma``). Set different from ``real_vol`` to study hedging at the
            wrong vol.
        real_vol: volatility of the simulated path (defaults to ``sigma``).
        return_samples: if True, also return the raw per-path P&L list.

    Returns a :class:`HedgeResult` (and the samples if requested).
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    hv = hedge_vol if hedge_vol is not None else sigma
    rv = real_vol if real_vol is not None else sigma

    rng = random.Random(seed)
    dt = t / n_steps
    drift = (b - 0.5 * rv * rv) * dt
    vol_step = rv * math.sqrt(dt)
    growth = math.exp(r * dt)   # cash growth per step

    premium = bsm_price(S, K, t, r, sigma, ot, b=b)

    samples = []
    for _ in range(n_paths):
        s = S
        tau = t
        # Short the option: collect premium as cash, buy delta shares.
        d = bsm_delta(s, K, tau, r, hv, ot, b=b)
        cash = premium - d * s
        shares = d

        for step in range(n_steps):
            z = rng.gauss(0.0, 1.0)
            s = s * math.exp(drift + vol_step * z)
            tau = t - (step + 1) * dt
            cash *= growth
            if tau > 1e-12:
                d_new = bsm_delta(s, K, tau, r, hv, ot, b=b)
            else:
                d_new = 0.0
            # Rebalance to d_new shares, funding the trade from cash.
            cash -= (d_new - shares) * s
            shares = d_new

        # At expiry: liquidate shares, pay the option payoff we are short.
        payoff = max(s - K, 0.0) if ot is OptionType.CALL else max(K - s, 0.0)
        pnl = cash + shares * s - payoff
        samples.append(pnl)

    n = len(samples)
    mean = sum(samples) / n
    var = sum((x - mean) ** 2 for x in samples) / (n - 1) if n > 1 else 0.0
    result = HedgeResult(mean_pnl=mean, std_pnl=math.sqrt(var),
                         min_pnl=min(samples), max_pnl=max(samples),
                         n_paths=n, n_steps=n_steps)
    if return_samples:
        return result, samples
    return result
