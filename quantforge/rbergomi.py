"""Rough Bergomi stochastic-volatility model (Bayer-Friz-Gatheral 2016).

The rough Bergomi model drives the instantaneous variance by a *rough* Gaussian
process -- a fractional Brownian motion with Hurst exponent ``H < 1/2`` -- which
reproduces the steep short-maturity volatility skew that classical diffusive
stochastic-vol models (Heston, SABR) miss. The spot and variance are

    dS_t / S_t = sqrt(V_t) dB_t,      B = rho W + sqrt(1 - rho^2) W_perp
    V_t = xi0 * exp(eta * Wtilde_t - 0.5 * eta^2 * t^{2H})

where ``Wtilde_t = sqrt(2H) int_0^t (t-s)^{H-1/2} dW_s`` is the Riemann-Liouville
fractional process (normalised so ``E[Wtilde_t^2] = t^{2H}``). The compensator
keeps the variance a martingale in expectation: ``E[V_t] = xi0`` for all ``t``.

Simulation uses the hybrid scheme of Bennedsen, Lunde & Pakkanen (2017) with
``kappa = 1``: the singular cell nearest each time point is sampled *exactly*
from the joint law of the Brownian increment and its kernel-weighted integral,
while the remaining cells use the optimally-placed discretised kernel. Pure
standard library; this is a validation-grade rough-vol smile generator, not a
production HFT pricer.
"""

import math
import random

from .montecarlo import MCResult, _summarize
from .bsm import OptionType, _coerce_type


def _hybrid_weights(H, dt, n):
    """Optimal discretised-kernel weights w_k for lags k = 2..n (BLP b*)."""
    alpha = H - 0.5
    w = [0.0, 0.0]  # k=0,1 unused (k=1 is the exact singular cell)
    if abs(alpha) < 1e-12:
        # H = 1/2: the kernel (t-s)^alpha is identically 1, so every far cell
        # contributes an unweighted Brownian increment (standard BM recovered).
        return w + [1.0] * (n - 1)
    for k in range(2, n + 1):
        b_star = ((k ** (alpha + 1.0) - (k - 1) ** (alpha + 1.0))
                  / (alpha + 1.0)) ** (1.0 / alpha)
        w.append((b_star * dt) ** alpha)
    return w


def rbergomi_paths(S, t, xi0, eta, H, rho, r=0.0, n_steps=100, n_paths=20_000,
                   antithetic=True, seed=None):
    """Yield discounted terminal spots under rough Bergomi (internal helper).

    Returns a list of terminal spot prices (already drifted by ``r``); the
    caller applies the payoff and discounting.
    """
    if not (0.0 < H < 1.0):
        raise ValueError("H must be in (0, 1)")
    if xi0 <= 0 or eta < 0:
        raise ValueError("xi0 must be positive and eta non-negative")
    if not (-1.0 <= rho <= 1.0):
        raise ValueError("rho must be in [-1, 1]")
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")

    dt = t / n_steps
    sqrt_dt = math.sqrt(dt)
    root2H = math.sqrt(2.0 * H)
    alpha = H - 0.5
    w = _hybrid_weights(H, dt, n_steps)

    # Singular-cell bivariate law of (increment dW, kernel integral Y):
    #   Var(dW)=dt, Cov(dW,Y)=dt^{a+1}/(a+1), Var(Y)=dt^{2a+1}/(2a+1).
    cov = dt ** (alpha + 1.0) / (alpha + 1.0)
    var_y = dt ** (2.0 * alpha + 1.0) / (2.0 * alpha + 1.0)
    y_beta = cov / dt                       # regression coeff of Y on dW
    y_resid = math.sqrt(max(var_y - cov * cov / dt, 0.0))
    rho_perp = math.sqrt(max(1.0 - rho * rho, 0.0))
    t_pow = [(i * dt) ** (2.0 * H) for i in range(n_steps + 1)]

    rng = random.Random(seed)
    terminals = []

    def one_path(z1s, z2s, zps):
        dW = [sqrt_dt * z for z in z1s]                 # Brownian increments
        Y = [y_beta * dW[i] + y_resid * z2s[i] for i in range(n_steps)]
        logS = math.log(S)
        v_prev = xi0                                    # V at t_0 (Wtilde_0 = 0)
        for i in range(n_steps):
            # Spot step over cell i uses the left-point variance (predictable).
            drift = (r - 0.5 * v_prev) * dt
            vol = math.sqrt(v_prev)
            logS += drift + vol * (rho * dW[i] + rho_perp * sqrt_dt * zps[i])
            # Volterra process at t_{i+1}: exact singular cell + far cells.
            wt = Y[i]
            for k in range(2, i + 2):
                wt += w[k] * dW[i - k + 1]
            wtilde = root2H * wt
            v_prev = xi0 * math.exp(eta * wtilde - 0.5 * eta * eta * t_pow[i + 1])
        return math.exp(logS)

    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        z1s = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        z2s = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        zps = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        terminals.append(one_path(z1s, z2s, zps))
        if antithetic:
            terminals.append(one_path([-z for z in z1s],
                                      [-z for z in z2s],
                                      [-z for z in zps]))
    return terminals


def rbergomi_price(S, K, t, xi0, eta, H, rho, r=0.0,
                   option_type=OptionType.CALL, n_steps=100, n_paths=20_000,
                   antithetic=True, seed=None) -> MCResult:
    """Monte Carlo a European option under rough Bergomi.

    Args:
        xi0: the (flat) forward variance curve level; ``sqrt(xi0)`` is the
            baseline vol.
        eta: vol-of-vol of the rough driver.
        H: Hurst exponent in (0, 1); ``H < 0.5`` is rough (steep short skew),
            ``H = 0.5`` recovers a standard lognormal-vol diffusion.
        rho: spot/vol correlation (negative gives the equity down-skew).

    Returns an :class:`MCResult`.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    sign = 1.0 if ot is OptionType.CALL else -1.0
    disc = math.exp(-r * t)
    terminals = rbergomi_paths(S, t, xi0, eta, H, rho, r, n_steps, n_paths,
                               antithetic, seed)
    samples = [disc * max(sign * (sT - K), 0.0) for sT in terminals]
    price, se = _summarize(samples)
    return MCResult(price=price, std_error=se, n_paths=len(samples))


def rbergomi_smile(S, strikes, t, xi0, eta, H, rho, r=0.0, n_steps=100,
                   n_paths=40_000, antithetic=True, seed=None):
    """Black-Scholes implied-vol smile a rough Bergomi model produces.

    Simulates one set of terminal spots and reprices every strike on it (common
    random numbers), then inverts each call price to its Black-Scholes implied
    vol. Returns ``(log_moneyness, vol)`` pairs sorted by strike on the forward
    ``F = S e^{r t}``. Rough dynamics (``H < 0.5``) with ``rho < 0`` give the
    steep negative short-maturity skew that motivates the model.
    """
    from .implied import implied_volatility

    F = S * math.exp(r * t)
    disc = math.exp(-r * t)
    terminals = rbergomi_paths(S, t, xi0, eta, H, rho, r, n_steps, n_paths,
                               antithetic, seed)
    out = []
    for K in sorted(strikes):
        samples = [disc * max(sT - K, 0.0) for sT in terminals]
        price, _ = _summarize(samples)
        try:
            iv = implied_volatility(price, S, K, t, r, OptionType.CALL, b=r)
        except ValueError:
            continue
        out.append((math.log(K / F), iv))
    return out
