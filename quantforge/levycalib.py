"""Calibrate exponential-Levy models to a market implied-vol smile.

Every model that plugs into the shared Carr-Madan engine (VG, NIG, Meixner,
CGMY) exposes a characteristic exponent ``psi(u; params)``. Given a set of
market Black-Scholes vols across strikes at one expiry, this module fits the
model parameters by least squares on implied vol, pricing each candidate smile
with a single Carr-Madan FFT strip (``carr_madan_strip``) rather than one
integral per strike, and optimising with the built-in Nelder-Mead.

A small registry maps a model name to its exponent, a smooth parameter
transform (so the optimiser works unconstrained while the raw parameters stay in
their valid region), and a seed. Pure standard library.
"""

import math

from .bsm import OptionType
from .optimize import nelder_mead
from .carrmadan import carr_madan_strip
from .implied import implied_volatility
from .variancegamma import _vg_psi
from .nig import _nig_psi
from .meixner import _meixner_psi
from .cgmy import _cgmy_psi


def _softplus(x):
    return math.log1p(math.exp(-abs(x))) + max(x, 0.0)


def _inv_softplus(y):
    y = max(y, 1e-8)
    return math.log(math.expm1(y)) if y < 30 else y


def _tanh_to(lo, hi, x):
    """Map R -> (lo, hi) smoothly."""
    return lo + 0.5 * (hi - lo) * (math.tanh(x) + 1.0)


def _inv_tanh_from(lo, hi, y):
    z = 2.0 * (y - lo) / (hi - lo) - 1.0
    z = max(min(z, 0.999), -0.999)
    return math.atanh(z)


# Each entry: (unpack(p)->raw params, pack(raw)->p, psi_builder(raw), seed_raw).
def _vg_spec():
    def unpack(p):
        return (_softplus(p[0]) + 1e-6, _softplus(p[1]) + 1e-6, p[2])
    def pack(raw):
        return [_inv_softplus(raw[0]), _inv_softplus(raw[1]), raw[2]]
    def build(raw):
        sigma, nu, theta = raw
        return lambda u: _vg_psi(u, sigma, nu, theta)
    return unpack, pack, build, (0.2, 0.3, -0.2)


def _nig_spec():
    def unpack(p):
        alpha = _softplus(p[0]) + 1e-3
        beta = _tanh_to(-alpha + 1e-6, alpha - 1e-6, p[1])
        delta = _softplus(p[2]) + 1e-6
        return (alpha, beta, delta)
    def pack(raw):
        alpha, beta, delta = raw
        return [_inv_softplus(alpha), _inv_tanh_from(-alpha, alpha, beta),
                _inv_softplus(delta)]
    def build(raw):
        alpha, beta, delta = raw
        return lambda u: _nig_psi(u, alpha, beta, delta)
    return unpack, pack, build, (15.0, -5.0, 0.5)


def _meixner_spec():
    def unpack(p):
        a = _softplus(p[0]) + 1e-3
        b = _tanh_to(-math.pi + 1e-4, math.pi - 1e-4, p[1])
        d = _softplus(p[2]) + 1e-6
        return (a, b, d)
    def pack(raw):
        a, b, d = raw
        return [_inv_softplus(a), _inv_tanh_from(-math.pi, math.pi, b),
                _inv_softplus(d)]
    def build(raw):
        a, b, d = raw
        return lambda u: _meixner_psi(u, a, b, d)
    return unpack, pack, build, (0.3, -0.3, 0.5)


def _cgmy_spec():
    def unpack(p):
        C = _softplus(p[0]) + 1e-6
        G = _softplus(p[1]) + 1e-3
        M = _softplus(p[2]) + 1e-3
        Y = _tanh_to(-2.0, 1.8, p[3])   # keep Y < 2, away from the pole
        return (C, G, M, Y)
    def pack(raw):
        C, G, M, Y = raw
        return [_inv_softplus(C), _inv_softplus(G), _inv_softplus(M),
                _inv_tanh_from(-2.0, 1.8, Y)]
    def build(raw):
        C, G, M, Y = raw
        return lambda u: _cgmy_psi(u, C, G, M, Y)
    return unpack, pack, build, (4.0, 5.0, 10.0, 0.5)


_SPECS = {
    "vg": _vg_spec, "nig": _nig_spec, "meixner": _meixner_spec, "cgmy": _cgmy_spec,
}


def levy_psi(model, params):
    """Return the characteristic exponent ``psi(u)`` for a named Levy model.

    ``model`` is one of ``"vg"``, ``"nig"``, ``"meixner"``, ``"cgmy"`` and
    ``params`` its raw parameter tuple (the same shape :func:`calibrate_levy_smile`
    returns). Handy for pricing or surface-building a model outside calibration.
    """
    model = model.lower()
    if model not in _SPECS:
        raise ValueError(f"unknown model {model!r}; choose from {sorted(_SPECS)}")
    _unpack, _pack, build, _seed = _SPECS[model]()
    return build(tuple(params))


def _strip_vols(S, t, r, q, psi, strikes, cm_alpha):
    """Model implied vols at ``strikes`` from one Carr-Madan FFT strip."""
    grid_K, grid_C = carr_madan_strip(S, t, r, q, psi, alpha=cm_alpha)
    n = len(grid_K)
    out = []
    for K in strikes:
        # Locate K in the (increasing) grid and linearly interpolate the call.
        lo, hi = 0, n - 1
        while lo < hi - 1:
            mid = (lo + hi) // 2
            if grid_K[mid] <= K:
                lo = mid
            else:
                hi = mid
        x0, x1 = grid_K[lo], grid_K[hi]
        y0, y1 = grid_C[lo], grid_C[hi]
        c = y0 + (y1 - y0) * (K - x0) / (x1 - x0) if x1 > x0 else y0
        try:
            iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=r - q)
        except ValueError:
            iv = float("nan")
        out.append(iv)
    return out


def calibrate_levy_smile(model, S, t, r, strikes, market_vols, q=0.0,
                         cm_alpha=1.5, max_iter=4000):
    """Fit a Levy model to a one-expiry market smile by least squares on vol.

    Args:
        model: one of ``"vg"``, ``"nig"``, ``"meixner"``, ``"cgmy"``.
        strikes, market_vols: matching sequences of strikes and Black-Scholes
            implied vols at expiry ``t``.
        cm_alpha: Carr-Madan damping used when pricing each candidate.

    Returns ``(params, rmse)`` where ``params`` is the fitted raw-parameter tuple
    for the chosen model and ``rmse`` is the root-mean-square implied-vol error.
    """
    model = model.lower()
    if model not in _SPECS:
        raise ValueError(f"unknown model {model!r}; choose from {sorted(_SPECS)}")
    strikes = [float(k) for k in strikes]
    market_vols = [float(v) for v in market_vols]
    n = len(strikes)
    if n != len(market_vols) or n < 3:
        raise ValueError("need at least three matching (strike, vol) points")

    unpack, pack, build, seed = _SPECS[model]()
    x0 = pack(seed)

    def objective(p):
        raw = unpack(p)
        try:
            model_vols = _strip_vols(S, t, r, q, build(raw), strikes, cm_alpha)
        except (ValueError, OverflowError, ZeroDivisionError):
            return 1e6
        err = 0.0
        for mv, iv in zip(market_vols, model_vols):
            if not math.isfinite(iv):
                return 1e6
            err += (iv - mv) ** 2
        return err

    best_p, best_f = nelder_mead(objective, x0, step=0.4, max_iter=max_iter,
                                 tol=1e-14)
    return unpack(best_p), math.sqrt(best_f / n)
