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
from .carrmadan import _fft


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


def _next_pow2(n):
    p = 1
    while p < n:
        p <<= 1
    return p


def _build_far_kernel(w, n_steps):
    """Precompute the FFT of the hybrid far-cell weight kernel (cached per run).

    The far-cell Volterra contribution ``F[i] = sum_j dW[j] g[i-j]`` with
    ``g[m] = w[m+1]`` is a convolution; returns ``(m_fft, G_fft)`` where
    ``G_fft`` is the kernel spectrum reused across every path.
    """
    m_fft = _next_pow2(2 * n_steps)
    g = [0.0] * m_fft
    for mm in range(1, n_steps):
        g[mm] = w[mm + 1]
    return m_fft, _fft(g)


def _far_sums_fft(dW, n_steps, m_fft, G_fft):
    """Far-cell sums F[i] for i=0..n_steps-1 via cached-kernel FFT convolution."""
    padded = list(dW) + [0.0] * (m_fft - n_steps)
    DW_fft = _fft(padded)
    prod = [DW_fft[i] * G_fft[i] for i in range(m_fft)]
    conv = _fft(prod, inverse=True)
    return [conv[i].real for i in range(n_steps)]


def rbergomi_paths(S, t, xi0, eta, H, rho, r=0.0, n_steps=100, n_paths=20_000,
                   antithetic=True, seed=None, fast="auto"):
    """Yield discounted terminal spots under rough Bergomi (internal helper).

    Returns a list of terminal spot prices (already drifted by ``r``); the
    caller applies the payoff and discounting.

    The Volterra convolution -- the ``O(n^2)`` inner sum of the hybrid-scheme
    far-cell weights against the Brownian increments -- can be evaluated with an
    FFT (``O(n log n)``): the weight kernel's FFT is precomputed once and reused
    across every path, so only one forward and one inverse transform of the
    increments happen per path. ``fast="auto"`` (default) uses the FFT when
    ``n_steps >= 200`` (where it beats the direct loop despite pure-Python FFT
    overhead) and the direct double loop otherwise; ``fast=True``/``False`` force
    the choice. The two paths agree to ~1e-12.
    """
    if fast == "auto":
        fast = n_steps >= 200
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

    # Kernel caching: the far-cell contribution to the Volterra process,
    #   F[i] = sum_{k=2}^{i+1} w[k] dW[i-k+1] = sum_j dW[j] g[i-j],  g[m]=w[m+1],
    # is a discrete convolution of the increments with a fixed weight kernel.
    # Precompute the kernel's FFT once so each path costs O(n log n), not O(n^2).
    if fast and n_steps > 1:
        m_fft, G_fft = _build_far_kernel(w, n_steps)
    else:
        m_fft = 0
        G_fft = None

    def one_path(z1s, z2s, zps):
        dW = [sqrt_dt * z for z in z1s]                 # Brownian increments
        Y = [y_beta * dW[i] + y_resid * z2s[i] for i in range(n_steps)]
        far = (_far_sums_fft(dW, n_steps, m_fft, G_fft)
               if (fast and n_steps > 1) else None)
        logS = math.log(S)
        v_prev = xi0                                    # V at t_0 (Wtilde_0 = 0)
        for i in range(n_steps):
            # Spot step over cell i uses the left-point variance (predictable).
            drift = (r - 0.5 * v_prev) * dt
            vol = math.sqrt(v_prev)
            logS += drift + vol * (rho * dW[i] + rho_perp * sqrt_dt * zps[i])
            # Volterra process at t_{i+1}: exact singular cell + far cells.
            if far is not None:
                wt = Y[i] + far[i]
            else:
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


def _rbergomi_w_stats(S, t, xi0, eta, H, rho, r, n_steps, n_paths,
                      antithetic, seed, fast="auto"):
    """Per-path W-measurable statistics for the conditional estimator.

    Conditioning on the volatility-driving Brownian motion ``W``, the terminal
    log-spot is Gaussian, so only two path functionals are needed to price any
    strike analytically:

        I1 = int_0^t sqrt(V_s) dW_s   (the rho-correlated stochastic integral)
        QV = int_0^t V_s ds           (the realised quadratic variation)

    Returns a list of ``(I1, QV)`` pairs (one per simulated W path).
    """
    if not (0.0 < H < 1.0):
        raise ValueError("H must be in (0, 1)")
    if xi0 <= 0 or eta < 0:
        raise ValueError("xi0 must be positive and eta non-negative")
    if not (-1.0 <= rho <= 1.0):
        raise ValueError("rho must be in [-1, 1]")
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    if fast == "auto":
        fast = n_steps >= 200

    dt = t / n_steps
    sqrt_dt = math.sqrt(dt)
    root2H = math.sqrt(2.0 * H)
    alpha = H - 0.5
    w = _hybrid_weights(H, dt, n_steps)
    cov = dt ** (alpha + 1.0) / (alpha + 1.0)
    var_y = dt ** (2.0 * alpha + 1.0) / (2.0 * alpha + 1.0)
    y_beta = cov / dt
    y_resid = math.sqrt(max(var_y - cov * cov / dt, 0.0))
    t_pow = [(i * dt) ** (2.0 * H) for i in range(n_steps + 1)]

    if fast and n_steps > 1:
        m_fft, G_fft = _build_far_kernel(w, n_steps)
    else:
        m_fft = 0
        G_fft = None

    rng = random.Random(seed)
    stats = []

    def one_path(z1s, z2s):
        dW = [sqrt_dt * z for z in z1s]
        Y = [y_beta * dW[i] + y_resid * z2s[i] for i in range(n_steps)]
        far = (_far_sums_fft(dW, n_steps, m_fft, G_fft)
               if (fast and n_steps > 1) else None)
        v_prev = xi0
        I1 = 0.0
        QV = 0.0
        for i in range(n_steps):
            vol = math.sqrt(v_prev)
            I1 += vol * dW[i]        # left-point sqrt(V) dW
            QV += v_prev * dt        # left-point V ds
            if far is not None:
                wt = Y[i] + far[i]
            else:
                wt = Y[i]
                for k in range(2, i + 2):
                    wt += w[k] * dW[i - k + 1]
            wtilde = root2H * wt
            v_prev = xi0 * math.exp(eta * wtilde - 0.5 * eta * eta * t_pow[i + 1])
        return (I1, QV)

    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        z1s = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        z2s = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        stats.append(one_path(z1s, z2s))
        if antithetic:
            stats.append(one_path([-z for z in z1s], [-z for z in z2s]))
    return stats


def _conditional_call(S, K, t, r, rho, I1, QV):
    """Black-Scholes call price conditional on one W path (see _rbergomi_w_stats).

    Given ``I1`` and ``QV``, log S_T | W is Gaussian, so the conditional call is
    a Black-Scholes price with an effective spot ``S exp(rho I1 - rho^2 QV / 2)``
    and an effective variance ``(1 - rho^2) QV``.
    """
    from .bsm import price as bs_price
    S_cond = S * math.exp(rho * I1 - 0.5 * rho * rho * QV)
    var_eff = (1.0 - rho * rho) * QV
    if var_eff <= 0.0:
        # rho = +-1: fully correlated, no residual noise -> discounted intrinsic.
        fwd = S_cond * math.exp(r * t)
        return math.exp(-r * t) * max(fwd - K, 0.0)
    sigma_eff = math.sqrt(var_eff / t)
    return bs_price(S_cond, K, t, r, sigma_eff, OptionType.CALL, b=r)


def rbergomi_price_cv(S, K, t, xi0, eta, H, rho, r=0.0, n_steps=100,
                      n_paths=20_000, antithetic=True, seed=None) -> MCResult:
    """Rough Bergomi European call by the conditional (turbocharged) estimator.

    Instead of simulating the orthogonal spot noise and averaging noisy payoffs
    (:func:`rbergomi_price`), this conditions on the volatility-driving Brownian
    motion and integrates the orthogonal noise out with a Black-Scholes formula
    (McCrickerd & Pakkanen, 2018). Every path contributes a smooth conditional
    price, so the Monte Carlo standard error drops sharply for the same paths --
    typically several-fold, and more as ``|rho|`` shrinks.

    Only calls are provided directly; puts follow from put-call parity on the
    forward ``S e^{r t}``.
    """
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    stats = _rbergomi_w_stats(S, t, xi0, eta, H, rho, r, n_steps, n_paths,
                              antithetic, seed)
    samples = [_conditional_call(S, K, t, r, rho, I1, QV) for I1, QV in stats]
    price, se = _summarize(samples)
    return MCResult(price=price, std_error=se, n_paths=len(samples))


def rbergomi_smile_cv(S, strikes, t, xi0, eta, H, rho, r=0.0, n_steps=100,
                      n_paths=40_000, antithetic=True, seed=None):
    """Rough Bergomi implied-vol smile via the conditional estimator.

    Like :func:`rbergomi_smile` but prices each strike with the low-variance
    conditional call on a shared set of W paths, then inverts to a Black-Scholes
    vol. Returns ``(log_moneyness, vol)`` pairs sorted by strike.
    """
    from .implied import implied_volatility

    F = S * math.exp(r * t)
    stats = _rbergomi_w_stats(S, t, xi0, eta, H, rho, r, n_steps, n_paths,
                              antithetic, seed)
    out = []
    for K in sorted(strikes):
        samples = [_conditional_call(S, K, t, r, rho, I1, QV) for I1, QV in stats]
        price, _ = _summarize(samples)
        try:
            iv = implied_volatility(price, S, K, t, r, OptionType.CALL, b=r)
        except ValueError:
            continue
        out.append((math.log(K / F), iv))
    return out


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
