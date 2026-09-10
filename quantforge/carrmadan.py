"""Shared Carr-Madan Fourier pricer for exponential-Levy models.

Any model whose log-price is ``ln S_T = ln S + (r - q + omega) t + L_t`` for a
Levy process ``L`` with a known characteristic exponent ``psi`` (so
``E[e^{i u L_t}] = e^{t psi(u)}``) can be priced by the Carr-Madan (1999)
transform: damp the call price by ``e^{alpha k}`` to make it integrable, then
recover it from a single Fourier integral of the risk-neutral characteristic
function. The martingale drift correction is ``omega = -psi(-i)``.

This module factors that machinery out so each Levy model only supplies its
``psi`` (see :mod:`quantforge.cgmy`, :mod:`quantforge.nig`); the Gauss-Legendre
integration nodes are shared with the Heston pricer. Pure standard library.
"""

import cmath
import math

from .bsm import OptionType, _coerce_type
from .heston import _GL_NODES, _GL_WEIGHTS


def carr_madan_call(S, K, t, r, q, psi, alpha=1.5, upper=200.0):
    """Carr-Madan European call price for a Levy model with exponent ``psi``.

    Args:
        psi: callable ``psi(u)`` returning the complex Levy characteristic
            exponent, ``E[e^{i u L_t}] = exp(t * psi(u))``.
        alpha: damping factor (> 0); needs ``E[S_T^{alpha+1}] < inf``, i.e.
            ``psi(-(alpha+1) i)`` finite.
        upper: truncation of the Fourier integral.

    Returns the undiscounted-to-priced European call value.
    """
    x0 = math.log(S)
    lnK = math.log(K)
    disc = math.exp(-r * t)
    omega = -psi(-1j)   # martingale correction: E[S_T] = S e^{(r-q)t}

    def char_logspot(u):
        drift = x0 + (r - q + omega) * t
        return cmath.exp(1j * u * drift + t * psi(u))

    half = 0.5 * upper
    total = 0.0
    for node, w in zip(_GL_NODES, _GL_WEIGHTS):
        nu = half * (node + 1.0)
        if nu <= 0:
            nu = 1e-8
        u = nu - (alpha + 1.0) * 1j
        phi = char_logspot(u)
        denom = alpha * alpha + alpha - nu * nu + 1j * (2.0 * alpha + 1.0) * nu
        rho = disc * phi / denom
        total += w * (cmath.exp(-1j * nu * lnK) * rho).real
    return math.exp(-alpha * lnK) * half * total / math.pi


def _fft(x, inverse=False):
    """In-place-style radix-2 Cooley-Tukey FFT of a complex list.

    ``len(x)`` must be a power of two. Returns a new list. Pure standard library;
    used by the Carr-Madan strip to price a whole log-strike grid in one pass.
    """
    n = len(x)
    if n & (n - 1) != 0:
        raise ValueError("FFT length must be a power of two")
    a = list(x)
    # Bit-reversal permutation.
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            a[i], a[j] = a[j], a[i]
    # Butterflies.
    length = 2
    sign = 1.0 if inverse else -1.0
    while length <= n:
        ang = sign * 2.0 * math.pi / length
        wlen = cmath.exp(1j * ang)
        for start in range(0, n, length):
            w = 1.0 + 0j
            half = length >> 1
            for k in range(half):
                u = a[start + k]
                v = a[start + k + half] * w
                a[start + k] = u + v
                a[start + k + half] = u - v
                w *= wlen
        length <<= 1
    if inverse:
        a = [v / n for v in a]
    return a


def carr_madan_strip(S, t, r, q, psi, alpha=1.5, n_fft=4096, eta=0.25):
    """Price a whole strip of European calls in one FFT (Carr-Madan 1999).

    Returns ``(strikes, calls)`` on a log-strike grid centred on the forward.
    The Carr-Madan damped-call transform is sampled at ``n_fft`` frequency points
    spaced ``eta`` apart, Simpson-weighted, and inverted with a single radix-2
    FFT -- so the entire smile costs one transform instead of one Gauss-Legendre
    integral per strike. The log-strike spacing is ``lambda = 2 pi / (n_fft eta)``.

    Args:
        psi: characteristic exponent ``psi(u)`` (as in :func:`levy_price`).
        alpha: damping factor (> 0); needs ``psi(-(alpha+1) i)`` finite.
        n_fft: FFT length (power of two).
        eta: frequency-grid spacing; smaller = finer strikes over a wider range.
    """
    if n_fft & (n_fft - 1) != 0:
        raise ValueError("n_fft must be a power of two")
    x0 = math.log(S)
    disc = math.exp(-r * t)
    omega = -psi(-1j)
    lam = 2.0 * math.pi / (n_fft * eta)          # log-strike spacing
    b = 0.5 * n_fft * lam                         # log-strike grid half-width
    k0 = x0 + (r - q) * t                          # centre grid on the forward

    def char_logspot(u):
        return cmath.exp(1j * u * (x0 + (r - q + omega) * t) + t * psi(u))

    vals = []
    for j in range(n_fft):
        nu = j * eta
        u = nu - (alpha + 1.0) * 1j
        denom = alpha * alpha + alpha - nu * nu + 1j * (2.0 * alpha + 1.0) * nu
        rho = disc * char_logspot(u) / denom
        # Simpson weights (3 + (-1)^{j+1} - delta_{j0}) / 3, with the grid phase.
        simpson = (3.0 + (-1.0) ** (j + 1) - (1.0 if j == 0 else 0.0)) / 3.0
        phase = cmath.exp(1j * nu * (b - k0))
        vals.append(eta * phase * rho * simpson)

    fft_out = _fft(vals, inverse=False)
    strikes, calls = [], []
    for m in range(n_fft):
        ku = -b + k0 + lam * m                     # log-strike
        call = (math.exp(-alpha * ku) / math.pi) * fft_out[m].real
        strikes.append(math.exp(ku))
        calls.append(call)
    return strikes, calls


def carr_madan_smile_strip(S, t, r, q, psi, k_lo=-0.5, k_hi=0.5,
                           alpha=1.5, n_fft=4096, eta=0.25):
    """Implied-vol smile over a log-moneyness window from one Carr-Madan FFT.

    Runs :func:`carr_madan_strip` once, keeps the grid strikes whose forward
    log-moneyness ``ln(K / F)`` lies in ``[k_lo, k_hi]``, and inverts each call
    to a Black-Scholes implied vol. Returns ``(log_moneyness, vol)`` pairs sorted
    by strike -- the whole smile from a single transform.
    """
    from .implied import implied_volatility

    strikes, calls = carr_madan_strip(S, t, r, q, psi, alpha=alpha,
                                      n_fft=n_fft, eta=eta)
    F = S * math.exp((r - q) * t)
    out = []
    for K, c in zip(strikes, calls):
        lm = math.log(K / F)
        if lm < k_lo or lm > k_hi or c <= 0.0:
            continue
        try:
            iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=r - q)
        except ValueError:
            continue
        out.append((lm, iv))
    return out


def levy_price(S, K, t, r, q, psi, option_type=OptionType.CALL,
               alpha=1.5, upper=200.0) -> float:
    """Price a European call/put for a Levy model via Carr-Madan + parity."""
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if t == 0:
        return max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)
    call = carr_madan_call(S, K, t, r, q, psi, alpha=alpha, upper=upper)
    if ot is OptionType.CALL:
        return call
    return call - S * math.exp(-q * t) + K * math.exp(-r * t)
