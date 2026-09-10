"""Heston Monte Carlo via Andersen's Quadratic-Exponential (QE) scheme.

Simulating the Heston stochastic-variance model naively (a plain Euler step on
the CIR variance) leaks negative variances and biases prices badly at the
maturities and vol-of-vol levels that matter. Andersen (2008), "Simple and
Efficient Simulation of the Heston Stochastic Volatility Model", fixes this with
two ideas:

  * **QE variance update.** Given the conditional mean ``m`` and variance ``s2``
    of the CIR process one step ahead, match the first two moments to one of two
    distributions depending on ``psi = s2 / m^2``. For ``psi <= psi_c`` (low
    vol-of-vol) use a shifted squared-Gaussian ``v' = a (b + Z)^2``; for
    ``psi > psi_c`` use an exponential with an atom at zero. Both stay
    non-negative by construction and match the exact CIR moments.

  * **Martingale-corrected log-asset step.** The log-spot is advanced with
    Andersen's ``K0..K4`` constants (his eq. 33), which fold the exact
    variance-integral drift and the leverage term ``rho`` into the step so the
    discounted spot stays a martingale.

Pure standard library; uses the same ``MCResult`` and variance-reduction
conventions as :mod:`quantforge.montecarlo`. This is the natural Monte Carlo
cross-check for the Fourier :func:`quantforge.heston_price`.
"""

import math
import random

from .montecarlo import MCResult, _summarize
from .bsm import OptionType, _coerce_type

PSI_C = 1.5  # switching threshold between the two QE branches (Andersen)


def heston_qe_mc(S, K, t, r, v0, kappa, theta, xi, rho,
                 option_type=OptionType.CALL, q=0.0, n_steps=100,
                 n_paths=50_000, antithetic=True, seed=None,
                 gamma1=0.5) -> MCResult:
    """Price a European option under Heston by Andersen's QE Monte Carlo.

    Args:
        v0, kappa, theta, xi, rho: Heston parameters (initial variance, mean
            reversion speed, long variance, vol-of-vol, spot/vol correlation).
        q: continuous dividend yield.
        n_steps: time steps (the QE variance update is exact in its moments, so
            the residual bias is only in the asset integral; ~50-200 is plenty).
        gamma1: weight on the left variance endpoint in the variance integral
            (``gamma1 = 0.5`` is the central discretisation; ``gamma2`` is set
            to ``1 - gamma1``).

    Returns an :class:`MCResult`. Puts come from simulating the same paths and
    taking the put payoff (parity holds path-by-path at the terminal spot).
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if v0 < 0 or theta < 0 or xi < 0:
        raise ValueError("variance parameters must be non-negative")
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")

    dt = t / n_steps
    sign = 1.0 if ot is OptionType.CALL else -1.0
    disc = math.exp(-r * t)

    ekt = math.exp(-kappa * dt)
    gamma2 = 1.0 - gamma1

    # Andersen's K-constants for the log-asset step (his eq. 33), with the
    # deterministic (r - q) dt drift kept explicit in K0d.
    K0 = -rho * kappa * theta * dt / xi
    K1 = gamma1 * dt * (kappa * rho / xi - 0.5) - rho / xi
    K2 = gamma2 * dt * (kappa * rho / xi - 0.5) + rho / xi
    K3 = gamma1 * dt * (1.0 - rho * rho)
    K4 = gamma2 * dt * (1.0 - rho * rho)
    K0d = K0 + (r - q) * dt

    rng = random.Random(seed)

    def next_var(v, uz):
        """One QE variance step. uz is a fresh U(0,1); returns v_{t+dt} >= 0."""
        m = theta + (v - theta) * ekt
        s2 = (v * xi * xi * ekt / kappa) * (1.0 - ekt) \
            + (theta * xi * xi / (2.0 * kappa)) * (1.0 - ekt) ** 2
        if m <= 0.0:
            return 0.0
        psi = s2 / (m * m)
        if psi <= PSI_C:
            # Shifted squared-Gaussian branch.
            inv = 2.0 / psi
            b2 = inv - 1.0 + math.sqrt(inv) * math.sqrt(max(inv - 1.0, 0.0))
            a = m / (1.0 + b2)
            # Invert the uniform to a standard normal via the inverse CDF.
            zv = _norm_ppf(uz)
            b = math.sqrt(b2)
            return a * (b + zv) ** 2
        else:
            # Exponential-with-atom branch.
            p = (psi - 1.0) / (psi + 1.0)
            beta = (1.0 - p) / m
            if uz <= p:
                return 0.0
            return math.log((1.0 - p) / (1.0 - uz)) / beta

    def one_path(norms, unis):
        x = math.log(S)     # log-spot
        v = v0
        for z, u in zip(norms, unis):
            v_next = next_var(v, u)
            vol = math.sqrt(max(K3 * v + K4 * v_next, 0.0))
            x += K0d + K1 * v + K2 * v_next + vol * z
            v = v_next
        sT = math.exp(x)
        return disc * max(sign * (sT - K), 0.0)

    samples = []
    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        norms = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        unis = [rng.random() for _ in range(n_steps)]
        samples.append(one_path(norms, unis))
        if antithetic:
            # Flip the asset normals; reflect the variance uniforms about 0.5 so
            # the antithetic path uses the mirror-image variance draws too.
            samples.append(one_path([-z for z in norms],
                                     [1.0 - u for u in unis]))

    price, se = _summarize(samples)
    return MCResult(price=price, std_error=se, n_paths=len(samples))


def heston_cv_mc(S, K, t, r, v0, kappa, theta, xi, rho,
                 option_type=OptionType.CALL, q=0.0, n_steps=100,
                 n_paths=50_000, antithetic=True, seed=None,
                 gamma1=0.5) -> MCResult:
    """Heston QE Monte Carlo with the underlying as a control variate.

    Andersen's QE step is martingale-corrected, so the discounted terminal spot
    ``Y = e^{-r t} S_T`` has the *known* mean ``E[Y] = S0 e^{-q t}`` (the
    discounted forward). ``Y`` is strongly correlated with the option payoff, so
    the controlled estimator ``X - beta (Y - E[Y])`` with the regression-optimal
    ``beta = Cov(X, Y) / Var(Y)`` sharply cuts the standard error at no bias.
    Everything else matches :func:`heston_qe_mc` (same QE variance step,
    ``K0..K4`` asset constants, antithetic draws).

    Cross-checks the Fourier :func:`quantforge.heston_price` and reports a
    standard error well below :func:`heston_qe_mc` at equal path count.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if v0 < 0 or theta < 0 or xi < 0:
        raise ValueError("variance parameters must be non-negative")
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")

    dt = t / n_steps
    sign = 1.0 if ot is OptionType.CALL else -1.0
    disc = math.exp(-r * t)
    ekt = math.exp(-kappa * dt)
    gamma2 = 1.0 - gamma1

    K0 = -rho * kappa * theta * dt / xi
    K1 = gamma1 * dt * (kappa * rho / xi - 0.5) - rho / xi
    K2 = gamma2 * dt * (kappa * rho / xi - 0.5) + rho / xi
    K3 = gamma1 * dt * (1.0 - rho * rho)
    K4 = gamma2 * dt * (1.0 - rho * rho)
    K0d = K0 + (r - q) * dt
    ey = S * math.exp(-q * t)                 # E[disc * S_T] under Heston

    rng = random.Random(seed)

    def next_var(v, uz):
        m = theta + (v - theta) * ekt
        s2 = (v * xi * xi * ekt / kappa) * (1.0 - ekt) \
            + (theta * xi * xi / (2.0 * kappa)) * (1.0 - ekt) ** 2
        if m <= 0.0:
            return 0.0
        psi = s2 / (m * m)
        if psi <= PSI_C:
            inv = 2.0 / psi
            b2 = inv - 1.0 + math.sqrt(inv) * math.sqrt(max(inv - 1.0, 0.0))
            a = m / (1.0 + b2)
            zv = _norm_ppf(uz)
            b = math.sqrt(b2)
            return a * (b + zv) ** 2
        p = (psi - 1.0) / (psi + 1.0)
        beta = (1.0 - p) / m
        if uz <= p:
            return 0.0
        return math.log((1.0 - p) / (1.0 - uz)) / beta

    def one_path(norms, unis):
        x = math.log(S)
        v = v0
        for z, u in zip(norms, unis):
            v_next = next_var(v, u)
            vol = math.sqrt(max(K3 * v + K4 * v_next, 0.0))
            x += K0d + K1 * v + K2 * v_next + vol * z
            v = v_next
        sT = math.exp(x)
        return disc * max(sign * (sT - K), 0.0), disc * sT

    xs, ys = [], []
    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        norms = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        unis = [rng.random() for _ in range(n_steps)]
        xp, yp = one_path(norms, unis)
        if antithetic:
            xm, ym = one_path([-z for z in norms], [1.0 - u for u in unis])
            xs.append(0.5 * (xp + xm))
            ys.append(0.5 * (yp + ym))
        else:
            xs.append(xp)
            ys.append(yp)

    m = len(xs)
    xbar = sum(xs) / m
    ybar = sum(ys) / m
    cov = sum((xx - xbar) * (yy - ybar) for xx, yy in zip(xs, ys))
    vary = sum((yy - ybar) ** 2 for yy in ys)
    beta = cov / vary if vary > 0.0 else 0.0

    controlled = [xx - beta * (yy - ey) for xx, yy in zip(xs, ys)]
    price, se = _summarize(controlled)
    return MCResult(price=price, std_error=se, n_paths=len(controlled))


def _norm_ppf(u):
    """Inverse standard-normal CDF (Acklam's rational approximation).

    Accurate to ~1e-9 over the open unit interval; used to turn the QE branch's
    uniform draw into the standard normal its squared-Gaussian form needs.
    """
    if u <= 0.0:
        return -1e10
    if u >= 1.0:
        return 1e10
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow = 0.02425
    phigh = 1.0 - plow
    if u < plow:
        ql = math.sqrt(-2.0 * math.log(u))
        return (((((c[0] * ql + c[1]) * ql + c[2]) * ql + c[3]) * ql + c[4]) * ql + c[5]) \
            / ((((d[0] * ql + d[1]) * ql + d[2]) * ql + d[3]) * ql + 1.0)
    if u > phigh:
        ql = math.sqrt(-2.0 * math.log(1.0 - u))
        return -(((((c[0] * ql + c[1]) * ql + c[2]) * ql + c[3]) * ql + c[4]) * ql + c[5]) \
            / ((((d[0] * ql + d[1]) * ql + d[2]) * ql + d[3]) * ql + 1.0)
    ql = u - 0.5
    rr = ql * ql
    return (((((a[0] * rr + a[1]) * rr + a[2]) * rr + a[3]) * rr + a[4]) * rr + a[5]) * ql \
        / (((((b[0] * rr + b[1]) * rr + b[2]) * rr + b[3]) * rr + b[4]) * rr + 1.0)


def heston_pathwise_delta(S, K, t, r, v0, kappa, theta, xi, rho,
                          option_type=OptionType.CALL, q=0.0, n_steps=100,
                          n_paths=50_000, antithetic=True, seed=None) -> MCResult:
    """Heston delta by the pathwise method (spot enters multiplicatively).

    In the QE simulation the initial spot appears only as the additive constant
    ``ln S0`` in the terminal log-price, so ``S_T = S0 e^Y`` with ``Y``
    independent of ``S0``. The pathwise delta is therefore exact and simple:

        delta = e^{-r t} E[ 1_{S_T > K} S_T / S0 ]   (put: -1_{S_T < K}),

    reusing the same Andersen-QE variance/asset scheme as
    :func:`heston_qe_mc`. Cross-checks a finite-difference bump of the price.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    dt = t / n_steps
    ekt = math.exp(-kappa * dt)
    gamma1 = 0.5
    gamma2 = 0.5
    K0 = -rho * kappa * theta * dt / xi
    K1 = gamma1 * dt * (kappa * rho / xi - 0.5) - rho / xi
    K2 = gamma2 * dt * (kappa * rho / xi - 0.5) + rho / xi
    K3 = gamma1 * dt * (1.0 - rho * rho)
    K4 = gamma2 * dt * (1.0 - rho * rho)
    K0d = K0 + (r - q) * dt
    disc = math.exp(-r * t)
    call = ot is OptionType.CALL
    rng = random.Random(seed)

    def next_var(v, uz):
        m = theta + (v - theta) * ekt
        s2 = (v * xi * xi * ekt / kappa) * (1.0 - ekt) \
            + (theta * xi * xi / (2.0 * kappa)) * (1.0 - ekt) ** 2
        if m <= 0.0:
            return 0.0
        psi = s2 / (m * m)
        if psi <= PSI_C:
            inv = 2.0 / psi
            b2 = inv - 1.0 + math.sqrt(inv) * math.sqrt(max(inv - 1.0, 0.0))
            a = m / (1.0 + b2)
            return a * (math.sqrt(b2) + _norm_ppf(uz)) ** 2
        p = (psi - 1.0) / (psi + 1.0)
        beta = (1.0 - p) / m
        return 0.0 if uz <= p else math.log((1.0 - p) / (1.0 - uz)) / beta

    def one(norms, unis):
        x = math.log(S)
        v = v0
        for z, u in zip(norms, unis):
            vn = next_var(v, u)
            vol = math.sqrt(max(K3 * v + K4 * vn, 0.0))
            x += K0d + K1 * v + K2 * vn + vol * z
            v = vn
        sT = math.exp(x)
        itm = (sT > K) if call else (sT < K)
        if not itm:
            return 0.0
        return disc * (sT / S) if call else -disc * (sT / S)

    samples = []
    n = n_paths // 2 if antithetic else n_paths
    for _ in range(n):
        norms = [rng.gauss(0.0, 1.0) for _ in range(n_steps)]
        unis = [rng.random() for _ in range(n_steps)]
        samples.append(one(norms, unis))
        if antithetic:
            samples.append(one([-z for z in norms], [1.0 - u for u in unis]))
    m, se = _summarize(samples)
    return MCResult(price=m, std_error=se, n_paths=len(samples))
