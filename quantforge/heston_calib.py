"""Calibrate the Heston model to an implied-vol surface.

Fits the five Heston parameters ``(v0, kappa, theta, xi, rho)`` by least squares
on Black implied vol over a set of ``(expiry, strike, market_vol)`` quotes,
pricing each candidate with the Fourier :func:`quantforge.heston_price` and
inverting to a vol. Nelder-Mead optimises under a smooth reparametrization
keeping ``v0, kappa, theta, xi > 0`` and ``rho in (-1, 1)``. Optionally penalises
Feller-condition violations (``2 kappa theta < xi^2``) so the fit prefers a
well-behaved variance process. Pure standard library.
"""

import math

from .bsm import OptionType
from .optimize import nelder_mead
from .heston import heston_price
from .implied import implied_volatility


def _softplus(x):
    return math.log1p(math.exp(-abs(x))) + max(x, 0.0)


def _inv_softplus(y):
    y = max(y, 1e-8)
    return math.log(math.expm1(y)) if y < 30 else y


def _unpack(p):
    v0 = _softplus(p[0]) + 1e-6
    kappa = _softplus(p[1]) + 1e-3
    theta = _softplus(p[2]) + 1e-6
    xi = _softplus(p[3]) + 1e-4
    rho = math.tanh(p[4])
    return (v0, kappa, theta, xi, rho)


def _pack(params):
    v0, kappa, theta, xi, rho = params
    return [_inv_softplus(v0), _inv_softplus(kappa), _inv_softplus(theta),
            _inv_softplus(xi), math.atanh(max(min(rho, 0.999), -0.999))]


def calibrate_heston(S, r, quotes, q=0.0, initial=None, feller_weight=0.0,
                     max_iter=6000):
    """Fit the five Heston parameters to an implied-vol surface.

    Args:
        quotes: iterable of ``(expiry, strike, market_vol)`` Black implied-vol
            points.
        initial: optional ``(v0, kappa, theta, xi, rho)`` seed; an ATM-variance
            based guess is used otherwise.
        feller_weight: if > 0, adds ``feller_weight * max(0, xi^2 - 2 kappa
            theta)^2`` to the objective, nudging the fit toward the Feller
            condition (a strictly positive variance process).

    Returns ``(params, rmse)`` -- the fitted 5-tuple and the root-mean-square
    implied-vol error over the quotes.
    """
    quotes = [(float(t), float(k), float(v)) for t, k, v in quotes]
    if len(quotes) < 4:
        raise ValueError("need at least four (expiry, strike, vol) quotes")

    if initial is None:
        atm = min(quotes, key=lambda z: abs(z[1] - S))[2]
        var = atm * atm
        initial = (var, 2.0, var, 0.4, -0.6)
    x0 = _pack(initial)

    def objective(p):
        v0, kappa, theta, xi, rho = _unpack(p)
        err = 0.0
        for t, K, mv in quotes:
            try:
                c = heston_price(S, K, t, r, v0, kappa, theta, xi, rho,
                                 OptionType.CALL, q=q)
                iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=r - q)
            except (ValueError, OverflowError, ZeroDivisionError):
                return 1e6
            err += (iv - mv) ** 2
        if feller_weight > 0.0:
            viol = max(0.0, xi * xi - 2.0 * kappa * theta)
            err += feller_weight * viol * viol
        return err

    best, f = nelder_mead(objective, x0, step=0.3, max_iter=max_iter, tol=1e-14)
    params = _unpack(best)
    # Report the pure vol RMSE (excluding any Feller penalty).
    sse = 0.0
    for t, K, mv in quotes:
        c = heston_price(S, K, t, r, *params, OptionType.CALL, q=q)
        iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=r - q)
        sse += (iv - mv) ** 2
    return params, math.sqrt(sse / len(quotes))
