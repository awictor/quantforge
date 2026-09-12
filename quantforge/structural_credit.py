"""Merton (1974) structural credit model.

The firm's assets ``V`` follow geometric Brownian motion; equity is a call option
on the assets struck at the face value of debt ``D`` maturing at ``T``, and the
firm defaults if ``V_T < D``. This module prices the equity, the risky debt and
its credit spread, and the risk-neutral default probability and distance to
default. Pure standard library on top of :mod:`quantforge.bsm`.
"""

import math

from .bsm import call_price
from .mathfns import norm_cdf


def _d1_d2(asset_value, debt_face, r, sigma, t):
    vsqrt = sigma * math.sqrt(t)
    d1 = (math.log(asset_value / debt_face) + (r + 0.5 * sigma * sigma) * t) / vsqrt
    return d1, d1 - vsqrt


def equity_value(asset_value, debt_face, r, sigma, t):
    """Equity as a call on the firm's assets struck at the debt face value.

    ``E = call(V, K=D, T)`` -- shareholders own the residual after repaying debt,
    a call on the assets. Increases with asset value and volatility.
    """
    if asset_value <= 0 or debt_face <= 0 or t <= 0 or sigma <= 0:
        raise ValueError("asset_value, debt_face, t, sigma must be positive")
    return call_price(asset_value, debt_face, t, r, sigma, b=r)


def equity_volatility(asset_value, debt_face, r, asset_vol, t):
    """Equity volatility implied by the asset volatility (Merton).

    From Ito's lemma on the equity call, ``sigma_E = (V/E) N(d1) sigma_V`` -- the
    equity is a levered claim, so its volatility exceeds the asset volatility by
    the delta-elasticity factor ``(V/E) N(d1)``. Rises as leverage rises.
    """
    if asset_value <= 0 or debt_face <= 0 or t <= 0 or asset_vol <= 0:
        raise ValueError("inputs must be positive")
    d1, _ = _d1_d2(asset_value, debt_face, r, asset_vol, t)
    E = call_price(asset_value, debt_face, t, r, asset_vol, b=r)
    if E <= 0:
        raise ValueError("equity value must be positive")
    return asset_value * norm_cdf(d1) * asset_vol / E


def solve_asset_value_and_vol(equity_value_obs, equity_vol_obs, debt_face, r, t,
                              tol=1e-10, max_iter=500):
    """Recover the unobservable asset value and volatility (KMV two-equation solve).

    Given the observed equity value and equity volatility, jointly solves the
    Merton system

        E = call(V, D, T),   sigma_E = (V/E) N(d1) sigma_V

    for ``(V, sigma_V)`` by fixed-point iteration: invert the equity-call for ``V``
    at the current ``sigma_V``, then update ``sigma_V`` from the equity-vol
    relation. Returns ``(asset_value, asset_vol)``. Round-trips with
    :func:`equity_value` and :func:`equity_volatility`.
    """
    if equity_value_obs <= 0 or equity_vol_obs <= 0 or debt_face <= 0 or t <= 0:
        raise ValueError("inputs must be positive")
    # Initial guesses: assets ~ equity + debt, asset vol ~ equity vol scaled down.
    V = equity_value_obs + debt_face * math.exp(-r * t)
    sigma_V = equity_vol_obs * equity_value_obs / V
    for _ in range(max_iter):
        # Invert E = call(V, D, T) for V at the current sigma_V (bisection).
        lo, hi = debt_face * 1e-6, equity_value_obs + debt_face * 10.0
        for _b in range(200):
            mid = 0.5 * (lo + hi)
            e_mid = call_price(mid, debt_face, t, r, sigma_V, b=r)
            if abs(e_mid - equity_value_obs) < tol:
                break
            if e_mid < equity_value_obs:
                lo = mid
            else:
                hi = mid
        V_new = mid
        d1, _ = _d1_d2(V_new, debt_face, r, sigma_V, t)
        sigma_V_new = equity_vol_obs * equity_value_obs / (V_new * norm_cdf(d1))
        if abs(V_new - V) < tol and abs(sigma_V_new - sigma_V) < tol:
            V, sigma_V = V_new, sigma_V_new
            break
        V, sigma_V = V_new, sigma_V_new
    return V, sigma_V


def risk_neutral_default_probability(asset_value, debt_face, r, sigma, t):
    """Risk-neutral probability of default ``P(V_T < D) = Phi(-d2)``.

    The chance the assets end below the debt face at maturity under the pricing
    measure. Rises with leverage (``D/V``), volatility, and horizon.
    """
    if asset_value <= 0 or debt_face <= 0 or t <= 0 or sigma <= 0:
        raise ValueError("asset_value, debt_face, t, sigma must be positive")
    _, d2 = _d1_d2(asset_value, debt_face, r, sigma, t)
    return norm_cdf(-d2)


def distance_to_default(asset_value, debt_face, r, sigma, t):
    """Distance to default: ``d2`` (standard deviations of asset drift above debt).

    ``d2 = (ln(V/D) + (r - sigma^2/2) T) / (sigma sqrt(T))``. The number of
    standard deviations the log assets must fall to hit the default barrier;
    higher is safer, and ``Phi(-DD)`` is the default probability.
    """
    if asset_value <= 0 or debt_face <= 0 or t <= 0 or sigma <= 0:
        raise ValueError("asset_value, debt_face, t, sigma must be positive")
    _, d2 = _d1_d2(asset_value, debt_face, r, sigma, t)
    return d2


def risky_debt_value(asset_value, debt_face, r, sigma, t):
    """Value of the risky debt: assets minus equity (``V - E``).

    By the accounting identity ``V = E + D_risky`` the debt is the firm value less
    the equity call. Below the risk-free discounted face; the gap is the credit
    risk.
    """
    return asset_value - equity_value(asset_value, debt_face, r, sigma, t)


def credit_spread(asset_value, debt_face, r, sigma, t):
    """Continuously-compounded credit spread of the risky debt over the risk-free.

    ``spread = -ln(D_risky / D e^{-r T}) / T`` -- the yield pickup of the risky
    debt over discounting the face at the risk-free rate. Non-negative, zero in the
    no-default limit, and rising with leverage and volatility.
    """
    if t <= 0:
        raise ValueError("t must be positive")
    debt = risky_debt_value(asset_value, debt_face, r, sigma, t)
    riskfree_debt = debt_face * math.exp(-r * t)
    if debt <= 0:
        raise ValueError("risky debt value must be positive")
    yield_risky = -math.log(debt / debt_face) / t
    yield_free = -math.log(riskfree_debt / debt_face) / t
    return yield_risky - yield_free
