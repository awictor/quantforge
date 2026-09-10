"""Quanto options: a foreign-asset option settled in domestic currency.

A quanto pays the foreign-asset option payoff, in the foreign currency's
numbers, but converted to domestic currency at a *fixed* exchange rate. The
quanto adjustment shows up as a shift to the asset's cost of carry by
``-rho * sigma_S * sigma_fx``, where ``rho`` is the correlation between the
asset and the domestic/foreign FX rate. Everything else is Black-Scholes with
domestic discounting.

The effective carry for the foreign asset under the domestic risk-neutral
measure is

    b_q = r_foreign - q_asset - rho * sigma_S * sigma_fx,

and the option is discounted at the domestic rate ``r_domestic``. Setting
``rho = 0`` (or ``sigma_fx = 0``) removes the quanto adjustment.
"""

import math

from .bsm import price as bsm_price, OptionType, _coerce_type, _validate


def quanto_option(S, K, t, r_domestic, r_foreign, sigma_asset, sigma_fx, rho,
                  q_asset=0.0, option_type=OptionType.CALL) -> float:
    """Price a quanto option (fixed-FX foreign-asset option in domestic terms).

    Args:
        S, K: foreign-asset spot and strike (in foreign-asset units).
        r_domestic: domestic risk-free rate (used for discounting).
        r_foreign: foreign risk-free rate.
        sigma_asset: volatility of the foreign asset.
        sigma_fx: volatility of the domestic/foreign FX rate.
        rho: correlation between the asset and the FX rate.
        q_asset: dividend yield on the foreign asset.

    The price is in domestic currency per unit of the fixed exchange rate
    (multiply by the agreed FX level for the cash amount).
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma_asset)
    if not (-1.0 <= rho <= 1.0):
        raise ValueError("rho must be in [-1, 1]")
    if sigma_fx < 0:
        raise ValueError("sigma_fx must be non-negative")

    # Quanto-adjusted cost of carry; discount at the domestic rate.
    b_q = r_foreign - q_asset - rho * sigma_asset * sigma_fx
    return bsm_price(S, K, t, r_domestic, sigma_asset, ot, b=b_q)
