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

from .bsm import (
    price as bsm_price, OptionType, _coerce_type, _validate,
    delta as bsm_delta, gamma as bsm_gamma,
)


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


def quanto_option_greeks(S, K, t, r_domestic, r_foreign, sigma_asset, sigma_fx,
                         rho, q_asset=0.0, option_type=OptionType.CALL):
    """Greeks of a quanto option.

    The quanto price is a Black-Scholes price on the foreign asset with the
    quanto-adjusted carry ``b_q = r_foreign - q_asset - rho sigma_asset
    sigma_fx``, discounted domestically. The spot enters only through that BSM
    price, so ``delta`` and ``gamma`` are the exact BSM Greeks at ``b_q`` (no
    finite difference). ``vega`` (dV/dsigma_asset -- which also moves ``b_q``),
    ``fx_vega`` (dV/dsigma_fx, the quanto's exposure to FX volatility), and
    ``corr_vega`` (dV/drho) are central finite differences of the closed form.
    Returns a dict with ``price``, ``delta``, ``gamma``, ``vega``, ``fx_vega``,
    ``corr_vega``.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma_asset)
    if not (-1.0 <= rho <= 1.0):
        raise ValueError("rho must be in [-1, 1]")
    if sigma_fx < 0:
        raise ValueError("sigma_fx must be non-negative")
    b_q = r_foreign - q_asset - rho * sigma_asset * sigma_fx
    price = bsm_price(S, K, t, r_domestic, sigma_asset, ot, b=b_q)
    delta = bsm_delta(S, K, t, r_domestic, sigma_asset, ot, b=b_q)
    gamma = bsm_gamma(S, K, t, r_domestic, sigma_asset, b=b_q)

    def px(sa=sigma_asset, sfx=sigma_fx, rr=rho):
        return quanto_option(S, K, t, r_domestic, r_foreign, sa, sfx, rr,
                             q_asset, ot)

    hv = 1e-4
    vega = (px(sa=sigma_asset + hv) - px(sa=sigma_asset - hv)) / (2.0 * hv)
    # sigma_fx has a floor of 0, so use a one-sided difference near the boundary.
    if sigma_fx >= hv:
        fx_vega = (px(sfx=sigma_fx + hv) - px(sfx=sigma_fx - hv)) / (2.0 * hv)
    else:
        fx_vega = (px(sfx=sigma_fx + hv) - price) / hv
    hr = 1e-5
    corr_vega = (px(rr=min(rho + hr, 1.0 - 1e-9))
                 - px(rr=max(rho - hr, -1.0 + 1e-9))) / (2.0 * hr)
    return {"price": price, "delta": delta, "gamma": gamma, "vega": vega,
            "fx_vega": fx_vega, "corr_vega": corr_vega}


def compo_option(S, K, t, r_domestic, r_foreign, sigma_asset, sigma_fx, rho,
                 q_asset=0.0, option_type=OptionType.CALL) -> float:
    """Price a composite (compo) FX option: a foreign asset valued in domestic terms.

    Unlike a quanto (fixed FX), a compo option converts the foreign asset to
    domestic currency at the *floating* exchange rate, so the payoff is on the
    domestic-currency asset value ``X = S * FX``. Its volatility combines the
    asset and FX vols with their correlation:

        sigma_compo = sqrt(sigma_asset^2 + sigma_fx^2 + 2 rho sigma_asset sigma_fx)

    Both ``S`` and ``K`` are quoted in domestic currency (K is the domestic
    strike on the converted asset). Carry and discounting use the domestic rate;
    the foreign rate enters as the asset's dividend-like yield ``q_asset``.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma_asset)
    if not (-1.0 <= rho <= 1.0):
        raise ValueError("rho must be in [-1, 1]")
    if sigma_fx < 0:
        raise ValueError("sigma_fx must be non-negative")

    sigma_compo = math.sqrt(sigma_asset * sigma_asset + sigma_fx * sigma_fx
                            + 2.0 * rho * sigma_asset * sigma_fx)
    # Domestic-currency asset carries at the domestic rate less the asset yield.
    b = r_domestic - q_asset
    return bsm_price(S, K, t, r_domestic, sigma_compo, ot, b=b)


def compo_option_greeks(S, K, t, r_domestic, r_foreign, sigma_asset, sigma_fx,
                        rho, q_asset=0.0, option_type=OptionType.CALL):
    """Greeks of a composite (compo) FX option.

    A compo option is a Black-Scholes price on the domestic-currency asset value
    with the combined volatility
    ``sigma_compo = sqrt(sigma_asset^2 + sigma_fx^2 + 2 rho sigma_asset sigma_fx)``
    and carry ``b = r_domestic - q_asset``. The spot enters only through the BSM
    price, so ``delta`` and ``gamma`` are exact BSM Greeks (no finite difference).
    ``vega`` (dV/dsigma_asset), ``fx_vega`` (dV/dsigma_fx), and ``corr_vega``
    (dV/drho) are central finite differences of the closed form; unlike a quanto,
    a compo is *long* FX volatility (positive ``fx_vega``). Returns a dict with
    ``price``, ``delta``, ``gamma``, ``vega``, ``fx_vega``, ``corr_vega``.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma_asset)
    if not (-1.0 <= rho <= 1.0):
        raise ValueError("rho must be in [-1, 1]")
    if sigma_fx < 0:
        raise ValueError("sigma_fx must be non-negative")
    sigma_compo = math.sqrt(sigma_asset * sigma_asset + sigma_fx * sigma_fx
                            + 2.0 * rho * sigma_asset * sigma_fx)
    b = r_domestic - q_asset
    price = bsm_price(S, K, t, r_domestic, sigma_compo, ot, b=b)
    delta = bsm_delta(S, K, t, r_domestic, sigma_compo, ot, b=b)
    gamma = bsm_gamma(S, K, t, r_domestic, sigma_compo, b=b)

    def px(sa=sigma_asset, sfx=sigma_fx, rr=rho):
        return compo_option(S, K, t, r_domestic, r_foreign, sa, sfx, rr,
                            q_asset, ot)

    hv = 1e-4
    vega = (px(sa=sigma_asset + hv) - px(sa=sigma_asset - hv)) / (2.0 * hv)
    # sigma_fx has a floor of 0, so use a one-sided difference near the boundary.
    if sigma_fx >= hv:
        fx_vega = (px(sfx=sigma_fx + hv) - px(sfx=sigma_fx - hv)) / (2.0 * hv)
    else:
        fx_vega = (px(sfx=sigma_fx + hv) - price) / hv
    hr = 1e-5
    corr_vega = (px(rr=min(rho + hr, 1.0 - 1e-9))
                 - px(rr=max(rho - hr, -1.0 + 1e-9))) / (2.0 * hr)
    return {"price": price, "delta": delta, "gamma": gamma, "vega": vega,
            "fx_vega": fx_vega, "corr_vega": corr_vega}
