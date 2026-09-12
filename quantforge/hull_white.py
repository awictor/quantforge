"""Hull-White (extended Vasicek) one-factor short-rate model.

Hull-White adds a time-dependent drift to Vasicek,

    dr = (theta(t) - a r) dt + sigma dW,

so the model reprices *any* initial zero curve exactly while keeping analytic bond
prices. The fitted bond price is

    P(t, T) = A(t, T) exp(-B(t, T) r_t),   B = (1 - e^{-a(T-t)}) / a,

with ``A`` chosen from the market discount factors and the forward curve so that
``P(0, T)`` equals the input. This module works from an initial discount-factor
curve (a callable ``P0(T)``) and returns fitted zero prices and the ``B`` factor.
Pure standard library.
"""

import math


def hw_B(a, tau):
    """Hull-White ``B(t, T) = (1 - e^{-a*tau}) / a`` for ``tau = T - t``."""
    if a == 0.0:
        return tau
    return (1.0 - math.exp(-a * tau)) / a


def hw_zero_from_curve(P0, r0, a, sigma, t, T, f0=None, eps=1e-5):
    """Hull-White bond price ``P(t, T)`` fitted to an initial curve ``P0``.

    Parameters
    ----------
    P0 : callable
        Initial discount factor ``P0(T)`` observed today (``P0(0) = 1``).
    r0 : float
        Current short rate; for consistency it should equal the initial instant
        forward ``f(0,0)``.
    a, sigma : float
        Mean reversion and volatility (``a`` may be 0 for the Ho-Lee limit).
    t, T : float
        Valuation and maturity times, ``0 <= t <= T``.
    f0 : callable, optional
        Initial instantaneous forward ``f(0, t)``; defaults to a finite-difference
        of ``-ln P0``.
    eps : float
        Step for the forward finite difference.

    Returns
    -------
    float
        The fitted ``P(t, T)``. At ``t = 0`` it reproduces ``P0(T)`` exactly.
    """
    if sigma < 0 or t < 0 or T < t:
        raise ValueError("require sigma >= 0 and 0 <= t <= T")
    if f0 is None:
        def f0(u):
            pu = P0(max(u, 0.0))
            pv = P0(u + eps)
            return -(math.log(pv) - math.log(pu)) / eps

    B = hw_B(a, T - t)
    # Analytic A(t,T) that matches the initial curve (Hull-White 1994).
    lnA = (math.log(P0(T) / P0(t)) + B * f0(t)
           - (sigma * sigma) / (4.0 * a) * B * B * (1.0 - math.exp(-2.0 * a * t))
           if a != 0.0 else
           math.log(P0(T) / P0(t)) + B * f0(t)
           - 0.5 * sigma * sigma * t * B * B)
    return math.exp(lnA - B * r0)


def hw_bond_option(P0, a, sigma, t_option, t_bond, strike, is_call=True):
    """European option on a zero-coupon bond under Hull-White (analytic).

    Prices an option expiring at ``t_option`` on a zero maturing at ``t_bond``,
    struck at ``strike``, using the initial discount curve ``P0``. The forward
    bond price is log-normal with volatility

        sigma_P = sigma * B(a, t_bond - t_option) * sqrt((1 - e^{-2 a t_option}) / (2 a)),

    giving a Black-style formula in the discount factors ``P0(t_bond)`` and
    ``P0(t_option)`` (Jamshidian / Hull-White). Returns the option value today.
    Requires ``0 < t_option < t_bond``.
    """
    import math as _m
    from .mathfns import norm_cdf as _N

    if sigma < 0 or not (0.0 < t_option < t_bond):
        raise ValueError("require sigma >= 0 and 0 < t_option < t_bond")
    pT = P0(t_bond)
    pS = P0(t_option)
    if a != 0.0:
        var = (sigma * sigma) * (1.0 - _m.exp(-2.0 * a * t_option)) / (2.0 * a)
    else:
        var = sigma * sigma * t_option
    sig_p = hw_B(a, t_bond - t_option) * _m.sqrt(var)
    if sig_p <= 0.0:
        # Degenerate: intrinsic on the forward.
        fwd = pT / pS
        payoff = max(fwd - strike, 0.0) if is_call else max(strike - fwd, 0.0)
        return pS * payoff
    h = _m.log(pT / (pS * strike)) / sig_p + 0.5 * sig_p
    if is_call:
        return pT * _N(h) - strike * pS * _N(h - sig_p)
    return strike * pS * _N(-h + sig_p) - pT * _N(-h)
