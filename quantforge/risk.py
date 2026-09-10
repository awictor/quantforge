"""Portfolio Value-at-Risk and Expected Shortfall.

Three complementary estimators over a :class:`~quantforge.portfolio.Book`:

  * ``parametric_var`` — delta-gamma (Cornish-Fisher) VaR assuming a normal
    shock to a single underlying; captures convexity via the book gamma;
  * ``historical_var`` — empirical VaR/ES from a set of realized return
    scenarios, revaluing the delta-gamma P&L on each;
  * ``montecarlo_var`` — VaR/ES by full repricing of every position under
    simulated shocks (no delta-gamma approximation).

VaR is reported as a positive number: the loss (in currency) not expected to be
exceeded at the given confidence over the horizon. Expected Shortfall (a.k.a.
CVaR) is the mean loss conditional on breaching the VaR.

The parametric and historical estimators assume the book is driven by one
underlying (all positions share the same spot), which is the common
single-name / single-index desk case. ``montecarlo_var`` reprices exactly and
extends naturally once multi-asset correlation is added.
"""

import math
import random
from dataclasses import dataclass

from .mathfns import norm_ppf, norm_pdf
from .bsm import greeks
from .portfolio import price_book, Contract


@dataclass(frozen=True)
class VaRResult:
    var: float            # positive loss number
    expected_shortfall: float
    confidence: float
    horizon_days: float
    method: str


def _z_and_tail(confidence):
    # Loss quantile lives in the lower tail of P&L.
    alpha = 1.0 - confidence
    z = norm_ppf(alpha)          # negative
    return alpha, z


def parametric_var(book, sigma_annual, spot, confidence=0.99, horizon_days=1.0,
                   trading_days=252):
    """Delta-gamma VaR/ES via a Cornish-Fisher expansion.

    Args:
        book: a priced :class:`Book` (net.delta, net.gamma in spot terms).
        sigma_annual: annualized volatility of the underlying's returns.
        spot: current underlying price (to turn return shocks into price moves).
        confidence: e.g. 0.99.
        horizon_days / trading_days: scale vol to the VaR horizon.

    P&L ~= delta * dS + 0.5 * gamma * dS^2, with dS = spot * r_h and r_h normal
    with std ``sigma_h``. The gamma term makes P&L non-normal; Cornish-Fisher
    corrects the quantile using the P&L skewness.
    """
    sigma_h = sigma_annual * math.sqrt(horizon_days / trading_days)
    ds_sd = spot * sigma_h

    d = book.net.delta
    g = book.net.gamma

    # Moments of dP = d*dS + 0.5*g*dS^2 where dS ~ N(0, ds_sd^2).
    # E[dS]=0, E[dS^2]=s2, E[dS^3]=0, E[dS^4]=3 s2^2.
    s2 = ds_sd * ds_sd
    mean = 0.5 * g * s2
    var_pnl = (d * d) * s2 + 0.5 * (g * g) * s2 * s2  # Var(a X + b X^2), X~N(0,s2)
    std = math.sqrt(var_pnl) if var_pnl > 0 else 0.0

    # Third central moment of dP (drives skew from gamma):
    #   E[(dP-mean)^3] = 3 d^2 g s2^2 + g^3 s2^3  (odd delta terms vanish)
    m3 = 3.0 * (d * d) * g * (s2 ** 2) + (g ** 3) * (s2 ** 3)
    skew = m3 / (std ** 3) if std > 0 else 0.0

    alpha, z = _z_and_tail(confidence)
    # Cornish-Fisher adjusted standard normal quantile.
    z_cf = z + (z * z - 1.0) * skew / 6.0
    pnl_quantile = mean + std * z_cf         # a loss => negative
    var = max(-pnl_quantile, 0.0)

    # ES under the normal-with-CF approximation: ES = mean + std * phi(z)/alpha
    # (magnitude); use base z for the density (standard ES formula).
    es_pnl = mean - std * (norm_pdf(z) / alpha)
    es = max(-es_pnl, 0.0)
    return VaRResult(var=var, expected_shortfall=es, confidence=confidence,
                     horizon_days=horizon_days, method="parametric-delta-gamma")


def historical_var(book, return_scenarios, spot, confidence=0.99, horizon_days=1.0):
    """Historical VaR/ES: apply each realized return to a delta-gamma P&L.

    Args:
        return_scenarios: iterable of simple returns r (e.g. daily), already at
            the desired horizon.
    """
    d = book.net.delta
    g = book.net.gamma
    pnls = []
    for r in return_scenarios:
        ds = spot * r
        pnls.append(d * ds + 0.5 * g * ds * ds)
    return _empirical_var(pnls, confidence, horizon_days, "historical-delta-gamma")


def _empirical_var(pnls, confidence, horizon_days, method):
    if not pnls:
        raise ValueError("no scenarios provided")
    losses = sorted(-p for p in pnls)  # ascending losses; tail at the top
    n = len(losses)
    # Index of the VaR quantile among losses.
    idx = min(n - 1, max(0, int(math.ceil(confidence * n)) - 1))
    var = max(losses[idx], 0.0)
    tail = [x for x in losses if x >= losses[idx]]
    es = max(sum(tail) / len(tail), 0.0) if tail else var
    return VaRResult(var=var, expected_shortfall=es, confidence=confidence,
                     horizon_days=horizon_days, method=method)


def montecarlo_var(contracts, sigma_annual, confidence=0.99, horizon_days=1.0,
                   trading_days=252, n_paths=20_000, seed=None):
    """Full-repricing VaR/ES: shock the spot, reprice every option, measure P&L.

    Unlike the parametric/historical estimators this makes no delta-gamma
    approximation — each position is repriced under the shocked spot with time
    advanced by the horizon.
    """
    contracts = [c.normalized() for c in contracts]
    base = price_book(contracts)
    base_mv = base.net.market_value

    sigma_h = sigma_annual * math.sqrt(horizon_days / trading_days)
    dt = horizon_days / trading_days
    rng = random.Random(seed)

    pnls = []
    for _ in range(n_paths):
        z = rng.gauss(0.0, 1.0)
        shocked = []
        for c in contracts:
            new_S = c.S * math.exp(-0.5 * sigma_h * sigma_h + sigma_h * z)
            new_t = max(c.t - dt, 1e-9)
            shocked.append(Contract(S=new_S, K=c.K, t=new_t, r=c.r, sigma=c.sigma,
                                    option_type=c.option_type, b=c.b, qty=c.qty,
                                    multiplier=c.multiplier, label=c.label))
        mv = price_book(shocked).net.market_value
        pnls.append(mv - base_mv)

    return _empirical_var(pnls, confidence, horizon_days, "montecarlo-full-reval")
