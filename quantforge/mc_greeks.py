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


def asian_pathwise_vega(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                        n_steps=50, n_paths=100_000, antithetic=True,
                        seed=None) -> MCResult:
    """Pathwise vega of a fixed-strike arithmetic-average Asian call/put.

    Differentiates the payoff along each path with respect to ``sigma``. With
    ``S_i = S0 exp(sum (b - sig^2/2) dt + sig sqrt(dt) Z_j)``, the pathwise
    sensitivity of each monitored spot is ``dS_i/dsig = S_i * (W_i - sig t_i)``
    where ``W_i = sqrt(dt) sum_{j<=i} Z_j`` is the accumulated Brownian motion,
    so the average's derivative is ``dA/dsig = mean_i dS_i/dsig`` and the payoff
    derivative is ``disc * 1_{A>K} * dA/dsig`` (put: ``-1_{A<K}``). Lower
    variance than a bump for this Lipschitz payoff; the kink at ``A = K`` is a
    measure-zero set.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    dt = t / n_steps
    sqdt = math.sqrt(dt)
    disc = math.exp(-r * t)
    call = ot is OptionType.CALL
    rng = random.Random(seed)
    samples = []

    def one(zs):
        s = S
        w = 0.0                      # accumulated Brownian motion
        avg = 0.0
        davg = 0.0                   # d(sum S_i)/dsig
        for i, z in enumerate(zs):
            s *= math.exp((b - 0.5 * sigma * sigma) * dt + sigma * sqdt * z)
            w += sqdt * z
            tau = (i + 1) * dt
            avg += s
            davg += s * (w - sigma * tau)
        avg /= n_steps
        davg /= n_steps
        if call:
            return disc * davg if avg > K else 0.0
        return -disc * davg if avg < K else 0.0

    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        zs = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        samples.append(one(zs))
        if antithetic:
            samples.append(one([-z for z in zs]))
    m, se = _summarize(samples)
    return MCResult(price=m, std_error=se, n_paths=len(samples))
