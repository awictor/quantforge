"""Liability-driven investing: funding ratios, duration hedging, surplus risk.

A pension or insurer holds assets against a stream of future liability cashflows.
This module values the liabilities off a discount rate, tracks the funding ratio
(assets over liabilities) and the surplus, measures the liability's interest-rate
duration, sizes the duration hedge, and gives the surplus-at-risk under a
funded-status volatility. Pure standard library.
"""

import math

from .mathfns import norm_ppf


def liability_pv(cashflows, discount_rate):
    """Present value of a liability stream ``[(t, amount), ...]``.

    Continuously-compounded discounting ``sum_i CF_i e^{-r t_i}``.
    """
    if discount_rate <= -1.0:
        raise ValueError("discount_rate must exceed -100%")
    return sum(cf * math.exp(-discount_rate * t) for t, cf in cashflows)


def funding_ratio(assets, liabilities):
    """Funding ratio ``assets / liabilities`` (above 1 = surplus)."""
    if liabilities <= 0:
        raise ValueError("liabilities must be positive")
    return assets / liabilities


def surplus(assets, liabilities):
    """Plan surplus (deficit if negative): ``assets - liabilities``."""
    return assets - liabilities


def liability_duration(cashflows, discount_rate):
    """Macaulay duration of the liability stream (years).

    PV-weighted average cashflow time ``sum t_i PV_i / sum PV_i``. The interest-
    rate sensitivity the asset portfolio must match to immunize the surplus.
    """
    pv = liability_pv(cashflows, discount_rate)
    if pv <= 0:
        raise ValueError("liability PV must be positive")
    return sum(t * cf * math.exp(-discount_rate * t) for t, cf in cashflows) / pv


def hedge_ratio(asset_duration, asset_value, liability_duration_,
                liability_value):
    """Fraction of the liability dollar-duration hedged by the assets.

    ``(asset_duration * asset_value) / (liability_duration * liability_value)`` --
    the ratio of asset to liability dollar duration (DV01). One means the surplus
    is immune to a parallel rate move; below one leaves residual liability
    interest-rate risk.
    """
    denom = liability_duration_ * liability_value
    if denom <= 0:
        raise ValueError("liability dollar duration must be positive")
    return (asset_duration * asset_value) / denom


def required_hedge_duration(asset_value, liability_duration_, liability_value):
    """Asset duration that fully immunizes the surplus (hedge ratio = 1).

    ``liability_duration * liability_value / asset_value`` -- the duration the
    asset portfolio must carry so its dollar duration matches the liability's.
    """
    if asset_value <= 0:
        raise ValueError("asset_value must be positive")
    return liability_duration_ * liability_value / asset_value


def surplus_at_risk(assets, liabilities, surplus_volatility, confidence=0.95,
                    horizon=1.0):
    """Surplus-at-risk: the worst surplus loss at a confidence over a horizon.

    ``z * surplus_volatility * sqrt(horizon) * liabilities`` where ``z =
    Phi^{-1}(confidence)`` and ``surplus_volatility`` is the funded-status (surplus/
    liabilities) return volatility. A one-sided downside measure (positive number =
    potential shortfall), analogous to VaR for the plan surplus.
    """
    if not (0.5 < confidence < 1.0):
        raise ValueError("confidence must be in (0.5, 1)")
    if surplus_volatility < 0 or horizon < 0:
        raise ValueError("volatility and horizon must be non-negative")
    if liabilities <= 0:
        raise ValueError("liabilities must be positive")
    z = norm_ppf(confidence)
    return z * surplus_volatility * math.sqrt(horizon) * liabilities
