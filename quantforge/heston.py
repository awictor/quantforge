"""Heston (1993) stochastic-volatility pricing via the characteristic function.

The Heston model gives the variance its own mean-reverting diffusion:

    dS = (r - q) S dt + sqrt(v) S dW1
    dv = kappa (theta - v) dt + xi sqrt(v) dW2,   d<W1, W2> = rho dt

Parameters:
    v0    : initial variance (sigma^2 at t=0)
    kappa : mean-reversion speed of the variance
    theta : long-run variance level
    xi    : vol-of-vol (volatility of the variance process)
    rho   : correlation between spot and variance shocks

European options price semi-analytically: the call is recovered from two
probabilities P1, P2 obtained by integrating the model's characteristic
function. We use the Albrecher et al. "little Heston trap" formulation of the
characteristic function, which keeps the complex logarithm on the principal
branch and avoids the discontinuities of Heston's original form. The single
Fourier integral is evaluated with a fixed Gauss-Legendre rule (pure Python,
no SciPy).
"""

import cmath
import math

from .bsm import OptionType, _coerce_type


# 40-point Gauss-Legendre nodes/weights on [-1, 1], generated once via the
# Newton method on the Legendre polynomial. Enough for smooth Heston integrands.
def _gauss_legendre(n):
    nodes = [0.0] * n
    weights = [0.0] * n
    m = (n + 1) // 2
    for i in range(m):
        # Initial guess for the i-th root.
        x = math.cos(math.pi * (i + 0.75) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, 0.0
            for k in range(n):
                p0, p1 = ((2 * k + 1) * x * p0 - k * p1) / (k + 1), p0
            # p0 = P_n(x), derivative dp:
            dp = n * (x * p0 - p1) / (x * x - 1.0)
            dx = p0 / dp
            x -= dx
            if abs(dx) < 1e-15:
                break
        nodes[i] = -x
        nodes[n - 1 - i] = x
        w = 2.0 / ((1.0 - x * x) * dp * dp)
        weights[i] = w
        weights[n - 1 - i] = w
    return nodes, weights


_GL_NODES, _GL_WEIGHTS = _gauss_legendre(64)


def _char_func(phi, S, K, t, r, q, v0, kappa, theta, xi, rho, j):
    """Heston characteristic-function integrand pieces (little Heston trap).

    Returns the complex value of exp(C + D*v0 + i*phi*x) / (i*phi) used in the
    probability integrals, where x = ln(S).
    """
    x = math.log(S)
    a = kappa * theta
    if j == 1:
        u = 0.5
        b = kappa - rho * xi
    else:
        u = -0.5
        b = kappa

    rxi = rho * xi
    d = cmath.sqrt((rxi * phi * 1j - b) ** 2 - xi * xi * (2.0 * u * phi * 1j - phi * phi))
    # "Little trap" g2 = (b - rxi*i*phi - d) / (b - rxi*i*phi + d)
    g = (b - rxi * phi * 1j - d) / (b - rxi * phi * 1j + d)

    exp_dt = cmath.exp(-d * t)
    C = ((r - q) * phi * 1j * t
         + (a / (xi * xi)) * ((b - rxi * phi * 1j - d) * t
                              - 2.0 * cmath.log((1.0 - g * exp_dt) / (1.0 - g))))
    D = ((b - rxi * phi * 1j - d) / (xi * xi)) * ((1.0 - exp_dt) / (1.0 - g * exp_dt))

    f = cmath.exp(C + D * v0 + 1j * phi * x)
    return (cmath.exp(-1j * phi * math.log(K)) * f / (1j * phi)).real


def _probability(S, K, t, r, q, v0, kappa, theta, xi, rho, j, upper=200.0):
    """P_j via Gauss-Legendre integration of the characteristic function."""
    half = 0.5 * upper
    total = 0.0
    for node, w in zip(_GL_NODES, _GL_WEIGHTS):
        phi = half * (node + 1.0)  # map [-1,1] -> [0, upper]
        if phi <= 0:
            phi = 1e-8
        total += w * _char_func(phi, S, K, t, r, q, v0, kappa, theta, xi, rho, j)
    integral = half * total
    return 0.5 + integral / math.pi


def heston_price(S, K, t, r, v0, kappa, theta, xi, rho,
                 option_type=OptionType.CALL, q=0.0, upper=200.0) -> float:
    """Price a European option under the Heston model.

    Args:
        v0, kappa, theta, xi, rho: Heston parameters (see module docstring).
        q: continuous dividend yield.
        upper: truncation of the Fourier integral (200 is ample for typical
            parameters; raise for very long maturities or large xi).

    Returns the option price. Puts are obtained from put-call parity.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if v0 < 0 or theta < 0 or xi < 0:
        raise ValueError("variance parameters must be non-negative")

    if t == 0:
        intrinsic = max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)
        return intrinsic

    P1 = _probability(S, K, t, r, q, v0, kappa, theta, xi, rho, 1, upper)
    P2 = _probability(S, K, t, r, q, v0, kappa, theta, xi, rho, 2, upper)

    call = S * math.exp(-q * t) * P1 - K * math.exp(-r * t) * P2
    if ot is OptionType.CALL:
        return call
    # Put-call parity: P = C - S e^{-qt} + K e^{-rt}.
    return call - S * math.exp(-q * t) + K * math.exp(-r * t)
