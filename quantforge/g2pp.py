"""Two-factor Gaussian G2++ short-rate model (Brigo-Mercurio).

G2++ adds a second correlated mean-reverting factor to the Gaussian short-rate
family, so it can produce the humped and decorrelated term-structure moves a
one-factor (Hull-White / single Cheyette) model cannot:

    r(t) = x(t) + y(t) + phi(t),
    dx = -a x dt + sigma dW1,   dy = -b y dt + eta dW2,   d<W1,W2> = rho dt,

with ``phi`` fitted to the initial curve. Zero-coupon bonds are exponential-
affine off the initial discount curve,

    P(t,T) = P(0,T)/P(0,t) * exp( -B(a,t,T) x - B(b,t,T) y
                                  + 0.5 (V(t,T) - V(0,T) + V(0,t)) ),

with ``B(z,t,T) = (1 - e^{-z (T-t)}) / z`` and the closed-form variance ``V``.
A European option on a zero-coupon bond has an exact Gaussian formula because
``ln P(t,T)`` is Gaussian; caplets follow by the bond-put identity. Pure
standard library. ``rho = 0`` and one vanishing factor recovers Hull-White.
"""

import math

from .mathfns import norm_cdf


def _B(z, tau):
    """The affine coefficient ``(1 - e^{-z tau}) / z`` (limit tau as z -> 0)."""
    if abs(z) < 1e-12:
        return tau
    return (1.0 - math.exp(-z * tau)) / z


def g2pp_V(a, b, sigma, eta, rho, t, T):
    """The G2++ variance term V(t,T) (Brigo-Mercurio eq. 4.10)."""
    tau = T - t
    Va = (sigma * sigma / (a * a)) * (
        tau + (2.0 / a) * math.exp(-a * tau) - (1.0 / (2.0 * a)) * math.exp(-2.0 * a * tau)
        - 1.5 / a)
    Vb = (eta * eta / (b * b)) * (
        tau + (2.0 / b) * math.exp(-b * tau) - (1.0 / (2.0 * b)) * math.exp(-2.0 * b * tau)
        - 1.5 / b)
    Vab = (2.0 * rho * sigma * eta / (a * b)) * (
        tau + (math.exp(-a * tau) - 1.0) / a + (math.exp(-b * tau) - 1.0) / b
        - (math.exp(-(a + b) * tau) - 1.0) / (a + b))
    return Va + Vb + Vab


def zero_bond(P0T, P0t, x, y, a, b, sigma, eta, rho, t, T):
    """G2++ zero-coupon bond ``P(t,T)`` given the factor state ``(x, y)``."""
    Ba = _B(a, T - t)
    Bb = _B(b, T - t)
    A = (0.5 * (g2pp_V(a, b, sigma, eta, rho, t, T)
                - g2pp_V(a, b, sigma, eta, rho, 0.0, T)
                + g2pp_V(a, b, sigma, eta, rho, 0.0, t)))
    return (P0T / P0t) * math.exp(-Ba * x - Bb * y + A)


def _bond_vol(a, b, sigma, eta, rho, t, T):
    """Std dev of ln P(t,T) seen from today, for the bond-option formula."""
    Ba = _B(a, T - t)
    Bb = _B(b, T - t)
    # Terminal variances of x(t), y(t) and their covariance.
    var_x = sigma * sigma * (1.0 - math.exp(-2.0 * a * t)) / (2.0 * a)
    var_y = eta * eta * (1.0 - math.exp(-2.0 * b * t)) / (2.0 * b)
    cov_xy = (rho * sigma * eta / (a + b)) * (1.0 - math.exp(-(a + b) * t))
    var = Ba * Ba * var_x + Bb * Bb * var_y + 2.0 * Ba * Bb * cov_xy
    return math.sqrt(max(var, 0.0))


def bond_option(P0S, P0T, a, b, sigma, eta, rho, expiry, maturity, strike,
                is_call=True):
    """European option on a zero-coupon bond under G2++ (exact).

    Option expires at ``expiry`` on a bond maturing at ``maturity``, struck at
    ``strike``. ``P0S = P(0, expiry)``, ``P0T = P(0, maturity)``. Since
    ``ln P(expiry, maturity)`` is Gaussian the price is a Black-style formula on
    the forward bond ``P0T / P0S`` with the G2++ bond volatility.
    """
    sig_p = _bond_vol(a, b, sigma, eta, rho, expiry, maturity)
    if sig_p < 1e-14:
        fwd = P0T / P0S
        payoff = max(fwd - strike, 0.0) if is_call else max(strike - fwd, 0.0)
        return P0S * payoff
    d1 = (math.log(P0T / (strike * P0S)) + 0.5 * sig_p * sig_p) / sig_p
    d2 = d1 - sig_p
    if is_call:
        return P0T * norm_cdf(d1) - strike * P0S * norm_cdf(d2)
    return strike * P0S * norm_cdf(-d2) - P0T * norm_cdf(-d1)


def caplet(P0_reset, P0_pay, a, b, sigma, eta, rho, reset, pay, strike,
           notional=1.0):
    """Caplet on ``[reset, pay]`` under G2++ via the bond-put identity."""
    tau = pay - reset
    K_bond = 1.0 / (1.0 + strike * tau)
    put = bond_option(P0_reset, P0_pay, a, b, sigma, eta, rho, reset, pay,
                      K_bond, is_call=False)
    return notional * (1.0 + strike * tau) * put
