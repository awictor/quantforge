"""Rough-Heston option pricing via the fractional Riccati equation.

The rough-Heston model (El Euch & Rosenbaum, 2019) is the rough-volatility limit
of a Hawkes-driven Heston model: the variance is driven by a fractional kernel
with Hurst exponent ``H = alpha - 1/2 in (0, 1/2]``, giving the steep short-dated
skew of rough volatility while retaining an (almost) affine structure. Its
characteristic function is

    E[e^{i a X_t}] = exp( theta*kappa * I^1 h(a, .)(t) + v0 * I^{1-alpha} h(a, .)(t) ),

where ``h`` solves the *fractional* Riccati equation

    D^alpha h = 1/2(-a^2 - i a) + kappa(i a rho nu - 1) h + (kappa nu)^2/2 h^2,

with ``alpha = H + 1/2``. There is no closed form, so we solve the Riccati by the
fractional Adams (predictor-corrector / product-rectangle) scheme on a time grid
and integrate the fractional integrals by the same quadrature, then price by the
Gil-Pelaez two-probability Fourier inversion reusing the shared Gauss-Legendre
nodes. ``H = 1/2`` (``alpha = 1``) recovers the classical Heston model. Pure
standard library.
"""

import cmath
import math

from .bsm import OptionType, _coerce_type
from .heston import _GL_NODES, _GL_WEIGHTS


def _frac_riccati(a, t, H, kappa, theta, nu, rho, v0, n_grid):
    """Solve the fractional Riccati for argument ``a``; return the CF exponent.

    Uses the fractional rectangle (product-integration) scheme of Diethelm for
    ``D^alpha h = F(h)`` with ``alpha = H + 1/2``, and accumulates the two
    fractional integrals that build the characteristic-function exponent.
    """
    alpha = H + 0.5
    dt = t / n_grid
    # F(h) = 0.5(-a^2 - i a) + kappa(i a rho nu - 1) h + 0.5 (kappa nu)^2 h^2.
    c0 = 0.5 * (-a * a - 1j * a)
    c1 = kappa * (1j * a * rho * nu - 1.0)
    c2 = 0.5 * (kappa * nu) ** 2

    def F(h):
        return c0 + c1 * h + c2 * h * h

    # Fractional Adams predictor-corrector (Diethelm-Ford-Freed). The predictor
    # uses the product-rectangle weights b_j; the corrector adds the
    # product-trapezoid weights a_j, giving O(dt^{1+alpha}) accuracy so the
    # alpha = 1 case converges to the exact Heston Riccati ODE.
    g_alpha1 = math.gamma(alpha + 1.0)
    g_alpha2 = math.gamma(alpha + 2.0)
    b_coef = dt ** alpha / g_alpha1
    a_coef = dt ** alpha / g_alpha2

    def a_weight(j, k):
        # product-trapezoid weight for node j in the corrector at step k.
        if j == 0:
            return (k - 1) ** (alpha + 1.0) - (k - 1 - alpha) * k ** alpha
        if j == k:
            return 1.0
        return ((k - j + 1) ** (alpha + 1.0) + (k - j - 1) ** (alpha + 1.0)
                - 2.0 * (k - j) ** (alpha + 1.0))

    h = [0.0 + 0.0j] * (n_grid + 1)   # h[0] = 0
    Fh = [0.0 + 0.0j] * (n_grid + 1)
    Fh[0] = F(h[0])
    for k in range(1, n_grid + 1):
        # Predictor: product-rectangle sum.
        pred = 0.0 + 0.0j
        for j in range(0, k):
            w = b_coef * ((k - j) ** alpha - (k - 1 - j) ** alpha)
            pred += w * Fh[j]
        Fp = F(pred)
        # Corrector: product-trapezoid sum + the new endpoint via the predictor.
        corr = a_coef * a_weight(k, k) * Fp
        for j in range(0, k):
            corr += a_coef * a_weight(j, k) * Fh[j]
        h[k] = corr
        Fh[k] = F(h[k])

    # CF exponent = theta*kappa * I^1 h + v0 * I^{1-alpha} h, both by quadrature.
    # I^1 h (ordinary integral) via the trapezoid rule.
    I1 = 0.0 + 0.0j
    for k in range(1, n_grid + 1):
        I1 += 0.5 * (h[k] + h[k - 1]) * dt
    # I^{1-alpha} h at t: fractional integral of order 1-alpha (= 1/2 - H).
    beta = 1.0 - alpha
    if beta <= 1e-12:
        # alpha = 1 (Heston): I^0 h = h(t).
        Ifrac = h[n_grid]
    else:
        gb = math.gamma(beta + 1.0)
        cb = dt ** beta / gb
        Ifrac = 0.0 + 0.0j
        for j in range(0, n_grid):
            w = cb * ((n_grid - j) ** beta - (n_grid - 1 - j) ** beta)
            Ifrac += w * h[j + 1]
    return theta * kappa * I1 + v0 * Ifrac


def _rh_cf(a, S, t, r, q, H, kappa, theta, nu, rho, v0, n_grid):
    """Rough-Heston characteristic function of ln S_T at complex argument ``a``."""
    x = math.log(S)
    drift = 1j * a * (x + (r - q) * t)
    expo = _frac_riccati(a, t, H, kappa, theta, nu, rho, v0, n_grid)
    return cmath.exp(drift + expo)


def _rh_probability(S, K, t, r, q, H, kappa, theta, nu, rho, v0, j,
                    n_grid, upper):
    lnK = math.log(K)
    fwd = _rh_cf(-1j, S, t, r, q, H, kappa, theta, nu, rho, v0, n_grid)
    half = 0.5 * upper
    total = 0.0
    for node, w in zip(_GL_NODES, _GL_WEIGHTS):
        phi = half * (node + 1.0)
        if phi <= 0:
            phi = 1e-8
        if j == 1:
            cf = _rh_cf(phi - 1j, S, t, r, q, H, kappa, theta, nu, rho, v0,
                       n_grid) / fwd
        else:
            cf = _rh_cf(phi, S, t, r, q, H, kappa, theta, nu, rho, v0, n_grid)
        total += w * (cmath.exp(-1j * phi * lnK) * cf / (1j * phi)).real
    return 0.5 + half * total / math.pi


def rough_heston_price(S, K, t, r, v0, kappa, theta, nu, rho, H=0.1,
                       option_type=OptionType.CALL, q=0.0, n_grid=200,
                       upper=120.0) -> float:
    """Price a European option under the rough-Heston model.

    Args:
        v0, kappa, theta: initial variance, mean-reversion speed, long variance.
        nu: volatility of variance *relative to kappa* (El Euch-Rosenbaum
            convention). At ``H = 0.5`` this reduces to the classical Heston
            model with vol-of-vol ``xi = kappa * nu``.
        rho: spot/variance correlation (negative for the equity skew).
        H: Hurst exponent in (0, 0.5]; ``H = 0.5`` recovers classical Heston.
        n_grid: fractional-Riccati time-grid resolution (more = more accurate,
            O(n_grid^2) work per Fourier node).

    Puts use put-call parity.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if not (0.0 < H <= 0.5):
        raise ValueError("H must be in (0, 0.5]")
    if v0 < 0 or theta < 0 or nu < 0:
        raise ValueError("variance parameters must be non-negative")
    if t == 0:
        return max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)

    P1 = _rh_probability(S, K, t, r, q, H, kappa, theta, nu, rho, v0, 1,
                         n_grid, upper)
    P2 = _rh_probability(S, K, t, r, q, H, kappa, theta, nu, rho, v0, 2,
                         n_grid, upper)
    call = S * math.exp(-q * t) * P1 - K * math.exp(-r * t) * P2
    if ot is OptionType.CALL:
        return call
    return call - S * math.exp(-q * t) + K * math.exp(-r * t)


def rough_heston_smile(S, strikes, t, r, v0, kappa, theta, nu, rho, H=0.1,
                       q=0.0, n_grid=200):
    """Black-Scholes implied-vol smile the rough-Heston model produces.

    Returns ``(log_moneyness, vol)`` pairs sorted by strike on the forward
    ``F = S e^{(r-q) t}``. Small ``H`` steepens the short-dated skew beyond what
    classical Heston can reach.
    """
    from .implied import implied_volatility

    F = S * math.exp((r - q) * t)
    out = []
    for K in sorted(strikes):
        c = rough_heston_price(S, K, t, r, v0, kappa, theta, nu, rho, H,
                               OptionType.CALL, q=q, n_grid=n_grid)
        try:
            iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=r - q)
        except ValueError:
            continue
        out.append((math.log(K / F), iv))
    return out
