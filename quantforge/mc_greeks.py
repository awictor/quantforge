"""Monte Carlo Greeks by the likelihood-ratio and pathwise methods.

Two ways to get Greeks from a simulation without bumping and re-pricing:

  * **Likelihood-ratio (Malliavin-flavoured) weights.** Differentiate the
    density, not the payoff, so the estimator is ``E[payoff * weight]`` for a
    Greek-specific weight. It works even for discontinuous payoffs (digitals),
    where the pathwise method fails. For a one-step Black-Scholes terminal
    ``S_T = S0 exp((b - sigma^2/2) t + sigma sqrt(t) Z)`` the weights are

        delta:  Z / (S0 sigma sqrt(t))
        vega:   (Z^2 - 1)/sigma - Z sqrt(t)
        gamma:  (Z^2 - Z sigma sqrt(t) - 1) / (S0^2 sigma^2 t)

  * **Pathwise** derivative. Differentiate the payoff along the path; lower
    variance where it applies (smooth payoffs). Pathwise delta of a call is
    ``e^{-r t} 1_{S_T > K} S_T / S0``.

Both use the same standard-library RNG and antithetic variates as
:mod:`quantforge.montecarlo`, and are cross-checked against the closed-form
Black-Scholes Greeks. Pure standard library.
"""

import math
import random

from .bsm import OptionType, _coerce_type, _validate
from .montecarlo import _summarize, MCResult


def lr_greeks(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
              n_paths=100_000, antithetic=True, seed=None):
    """European delta, gamma, vega by the likelihood-ratio method.

    Returns a dict with ``price``, ``delta``, ``gamma``, ``vega`` (each a Monte
    Carlo mean) plus their ``*_se`` standard errors. Works for the discontinuous
    digital payoff too (the LR weights do not touch the payoff).
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    disc = math.exp(-r * t)
    vsqrt = sigma * math.sqrt(t)
    sign = 1.0 if ot is OptionType.CALL else -1.0
    rng = random.Random(seed)

    price, delta, gamma, vega = [], [], [], []

    def acc(z):
        sT = S * math.exp((b - 0.5 * sigma * sigma) * t + vsqrt * z)
        pay = disc * max(sign * (sT - K), 0.0)
        price.append(pay)
        delta.append(pay * z / (S * vsqrt))
        gamma.append(pay * (z * z - z * vsqrt - 1.0) / (S * S * sigma * sigma * t))
        vega.append(pay * ((z * z - 1.0) / sigma - z * math.sqrt(t)))

    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        z = rng.gauss(0.0, 1.0)
        acc(z)
        if antithetic:
            acc(-z)

    out = {}
    for name, samples in (("price", price), ("delta", delta),
                          ("gamma", gamma), ("vega", vega)):
        m, se = _summarize(samples)
        out[name] = m
        out[name + "_se"] = se
    return out


def pathwise_delta(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                   n_paths=100_000, antithetic=True, seed=None) -> MCResult:
    """European delta by the pathwise method (smooth-payoff estimator).

    Delta of a call is ``e^{-r t} 1_{S_T > K} S_T / S0`` (put: ``-1_{S_T < K}``).
    Lower variance than the likelihood-ratio delta for these Lipschitz payoffs.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    disc = math.exp(-r * t)
    vsqrt = sigma * math.sqrt(t)
    call = ot is OptionType.CALL
    rng = random.Random(seed)
    samples = []

    def one(z):
        sT = S * math.exp((b - 0.5 * sigma * sigma) * t + vsqrt * z)
        if call:
            return disc * (sT / S) if sT > K else 0.0
        return -disc * (sT / S) if sT < K else 0.0

    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        z = rng.gauss(0.0, 1.0)
        samples.append(one(z))
        if antithetic:
            samples.append(one(-z))
    m, se = _summarize(samples)
    return MCResult(price=m, std_error=se, n_paths=len(samples))


def lr_digital_delta(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                     cash=1.0, n_paths=200_000, antithetic=True,
                     seed=None) -> MCResult:
    """Delta of a cash-or-nothing digital by the likelihood-ratio method.

    The digital payoff is discontinuous, so pathwise delta is ill-defined, but
    the LR estimator ``E[payoff * Z/(S0 sigma sqrt t)]`` is fine.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    disc = math.exp(-r * t)
    vsqrt = sigma * math.sqrt(t)
    call = ot is OptionType.CALL
    rng = random.Random(seed)
    samples = []

    def one(z):
        sT = S * math.exp((b - 0.5 * sigma * sigma) * t + vsqrt * z)
        itm = (sT > K) if call else (sT < K)
        pay = disc * cash if itm else 0.0
        return pay * z / (S * vsqrt)

    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        z = rng.gauss(0.0, 1.0)
        samples.append(one(z))
        if antithetic:
            samples.append(one(-z))
    m, se = _summarize(samples)
    return MCResult(price=m, std_error=se, n_paths=len(samples))
