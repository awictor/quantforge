"""Single-factor Cheyette (quasi-Gaussian) short-rate model.

The Cheyette model is a Markovian representation of the Heath-Jarrow-Morton
framework with a separable volatility. For a single factor with mean-reversion
speed ``kappa`` and short-rate volatility ``sigma(t)`` the short rate is

    r(t) = f(0, t) + x(t),
    dx = (y - kappa x) dt + sigma dW,   dy = (sigma^2 - 2 kappa y) dt,

where ``x`` is the deviation from the initial forward curve and ``y`` is the
(deterministic-given-vol) auxiliary state that makes the pair Markovian. Zero
-coupon bonds reconstitute analytically off the initial discount curve:

    P(t, T) = P(0,T)/P(0,t) * exp( -x G(t,T) - 0.5 y G(t,T)^2 ),
    G(t,T) = (1 - e^{-kappa (T - t)}) / kappa.

With a constant ``sigma`` this factor coincides with Hull-White, so ``x(t)`` is
Gaussian with closed-form mean and variance and bond options / caplets have
exact formulas -- which this module provides and cross-checks against. Pure
standard library.
"""

import math

from .mathfns import norm_cdf


def cheyette_G(kappa, tau):
    """The Cheyette/Hull-White G-function ``(1 - e^{-kappa tau}) / kappa``."""
    if abs(kappa) < 1e-12:
        return tau
    return (1.0 - math.exp(-kappa * tau)) / kappa


def _x_variance(kappa, sigma, t):
    """Var[x(t)] for constant sigma: sigma^2 (1 - e^{-2 kappa t}) / (2 kappa)."""
    if abs(kappa) < 1e-12:
        return sigma * sigma * t
    return sigma * sigma * (1.0 - math.exp(-2.0 * kappa * t)) / (2.0 * kappa)


def cheyette_y(kappa, sigma, t):
    """The auxiliary state ``y(t)`` for constant sigma (= Var[x(t)])."""
    return _x_variance(kappa, sigma, t)


def zero_bond(P0T, P0t, x, y, kappa, t, T):
    """Cheyette zero-coupon bond ``P(t, T)`` given the state ``(x, y)``.

    ``P0T = P(0, T)`` and ``P0t = P(0, t)`` are today's discount factors.
    """
    G = cheyette_G(kappa, T - t)
    return (P0T / P0t) * math.exp(-x * G - 0.5 * y * G * G)


def bond_option(P0S, P0T, kappa, sigma, expiry, maturity, strike,
                is_call=True):
    """Price a European option on a zero-coupon bond under Cheyette (const sigma).

    Option expires at ``expiry`` (= ``t``) on a bond maturing at ``maturity``
    (= ``T``), struck at ``strike``. ``P0S = P(0, expiry)`` and
    ``P0T = P(0, maturity)`` are today's discount factors. This is the exact
    Hull-White bond-option formula (Cheyette with constant sigma coincides with
    Hull-White), used as the analytic anchor for the model.
    """
    G = cheyette_G(kappa, maturity - expiry)
    var_x = _x_variance(kappa, sigma, expiry)
    sigma_p = G * math.sqrt(var_x)      # lognormal vol of P(t, T) at t
    if sigma_p < 1e-14:
        fwd = P0T / P0S
        payoff = max(fwd - strike, 0.0) if is_call else max(strike - fwd, 0.0)
        return P0S * payoff
    d1 = (math.log(P0T / (strike * P0S)) + 0.5 * sigma_p * sigma_p) / sigma_p
    d2 = d1 - sigma_p
    if is_call:
        return P0T * norm_cdf(d1) - strike * P0S * norm_cdf(d2)
    return strike * P0S * norm_cdf(-d2) - P0T * norm_cdf(-d1)


def caplet(P0_reset, P0_pay, kappa, sigma, reset, pay, strike, notional=1.0):
    """Price a caplet under Cheyette (constant sigma) on ``[reset, pay]``.

    A caplet paying ``tau (L - strike)^+`` at ``pay`` (with ``L`` the simple
    forward rate and ``tau = pay - reset``) equals ``notional (1 + strike tau)``
    put options on the zero-coupon bond ``P(reset, pay)`` struck at
    ``1 / (1 + strike tau)`` (the standard caplet<->bond-put identity).
    """
    tau = pay - reset
    K_bond = 1.0 / (1.0 + strike * tau)
    put = bond_option(P0_reset, P0_pay, kappa, sigma, reset, pay, K_bond,
                      is_call=False)
    return notional * (1.0 + strike * tau) * put
