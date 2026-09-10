"""Calibrate the double-Heston model to an implied-vol surface.

Double Heston has ten parameters -- two variance factors each with
``(v0, kappa, theta, xi, rho)``. Fitting all ten to a full surface is what makes
the model able to match both the level/term-structure and the skew that a single
factor cannot. This module fits them by least squares on Black implied vol over a
set of ``(expiry, strike, market_vol)`` quotes, pricing each candidate with the
Fourier :func:`quantforge.double_heston_price` and inverting to a vol, optimised
with the built-in Nelder-Mead under a smooth reparametrization that keeps every
parameter in its valid region (variances and vol-of-vols positive, correlations
in ``(-1, 1)``). Pure standard library.
"""

import math

from .bsm import OptionType
from .optimize import nelder_mead
from .double_heston import double_heston_price
from .implied import implied_volatility


def _softplus(x):
    return math.log1p(math.exp(-abs(x))) + max(x, 0.0)


def _inv_softplus(y):
    y = max(y, 1e-8)
    return math.log(math.expm1(y)) if y < 30 else y


def _unpack(p):
    """Map 10 unconstrained reals to valid double-Heston parameters."""
    v01 = _softplus(p[0]) + 1e-6
    kappa1 = _softplus(p[1]) + 1e-3
    theta1 = _softplus(p[2]) + 1e-6
    xi1 = _softplus(p[3]) + 1e-4
    rho1 = math.tanh(p[4])
    v02 = _softplus(p[5]) + 1e-6
    kappa2 = _softplus(p[6]) + 1e-3
    theta2 = _softplus(p[7]) + 1e-6
    xi2 = _softplus(p[8]) + 1e-4
    rho2 = math.tanh(p[9])
    return (v01, kappa1, theta1, xi1, rho1, v02, kappa2, theta2, xi2, rho2)


def _pack(params):
    v01, k1, th1, xi1, rho1, v02, k2, th2, xi2, rho2 = params
    return [_inv_softplus(v01), _inv_softplus(k1), _inv_softplus(th1),
            _inv_softplus(xi1), math.atanh(max(min(rho1, 0.999), -0.999)),
            _inv_softplus(v02), _inv_softplus(k2), _inv_softplus(th2),
            _inv_softplus(xi2), math.atanh(max(min(rho2, 0.999), -0.999))]


def calibrate_double_heston(S, r, quotes, q=0.0, initial=None, max_iter=6000):
    """Fit the ten double-Heston parameters to an implied-vol surface.

    Args:
        quotes: iterable of ``(expiry, strike, market_vol)`` Black implied-vol
            points.
        initial: optional starting 10-tuple ``(v01, kappa1, theta1, xi1, rho1,
            v02, kappa2, theta2, xi2, rho2)``; a sensible two-scale seed is used
            otherwise (a fast- and a slow-reverting factor).

    Returns ``(params, rmse)`` -- the fitted 10-tuple and the root-mean-square
    implied-vol error over the quotes.
    """
    quotes = [(float(t), float(k), float(v)) for t, k, v in quotes]
    if len(quotes) < 6:
        raise ValueError("need at least six (expiry, strike, vol) quotes")

    if initial is None:
        # Split ATM variance across a fast and a slow factor.
        atm = min(quotes, key=lambda z: abs(z[1] - S))[2]
        var = atm * atm
        initial = (0.5 * var, 3.0, 0.5 * var, 0.3, -0.6,
                   0.5 * var, 0.5, 0.5 * var, 0.2, -0.5)
    x0 = _pack(initial)

    def objective(p):
        pr = _unpack(p)
        err = 0.0
        for t, K, mv in quotes:
            try:
                c = double_heston_price(S, K, t, r, *pr, OptionType.CALL, q=q)
                iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=r - q)
            except (ValueError, OverflowError, ZeroDivisionError):
                return 1e6
            err += (iv - mv) ** 2
        return err

    best, f = nelder_mead(objective, x0, step=0.3, max_iter=max_iter, tol=1e-14)
    return _unpack(best), math.sqrt(f / len(quotes))
