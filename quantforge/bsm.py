"""Black-Scholes-Merton pricing and analytic Greeks.

The generalized BSM model prices European options on an asset with a
continuous cost-of-carry rate ``b``:

    b = r          -> Black-Scholes (non-dividend stock)
    b = r - q      -> Merton (continuous dividend yield q)
    b = 0          -> Black-76 (option on a future, discounted at r)
    b = r - r_f    -> Garman-Kohlhagen (FX, foreign rate r_f)

All functions take the carry ``b`` as a keyword defaulting to ``r`` so the
common stock case needs no extra arguments.
"""

import math
from dataclasses import dataclass
from enum import Enum

from .mathfns import norm_cdf, norm_pdf


class OptionType(str, Enum):
    CALL = "call"
    PUT = "put"


def _coerce_type(option_type) -> OptionType:
    if isinstance(option_type, OptionType):
        return option_type
    t = str(option_type).strip().lower()
    if t in ("c", "call"):
        return OptionType.CALL
    if t in ("p", "put"):
        return OptionType.PUT
    raise ValueError(f"unknown option type: {option_type!r}")


def _validate(S, K, t, sigma):
    if S <= 0:
        raise ValueError("spot S must be positive")
    if K <= 0:
        raise ValueError("strike K must be positive")
    if t < 0:
        raise ValueError("time to expiry t must be non-negative")
    if sigma < 0:
        raise ValueError("volatility sigma must be non-negative")


def _d1_d2(S, K, t, r, sigma, b):
    vsqrt = sigma * math.sqrt(t)
    d1 = (math.log(S / K) + (b + 0.5 * sigma * sigma) * t) / vsqrt
    d2 = d1 - vsqrt
    return d1, d2


def price(S, K, t, r, sigma, option_type=OptionType.CALL, b=None) -> float:
    """Price a European option under the generalized BSM model.

    Args:
        S: spot price of the underlying.
        K: strike price.
        t: time to expiry in years.
        r: continuously-compounded risk-free rate.
        sigma: annualized volatility.
        option_type: CALL or PUT (also accepts "call"/"put"/"c"/"p").
        b: cost of carry. Defaults to ``r`` (non-dividend stock).
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r

    # Handle degenerate limits cleanly (expiry or zero-vol -> intrinsic).
    if t == 0 or sigma == 0:
        fwd = S * math.exp(b * t)
        disc = math.exp(-r * t)
        if ot is OptionType.CALL:
            return disc * max(fwd - K, 0.0)
        return disc * max(K - fwd, 0.0)

    d1, d2 = _d1_d2(S, K, t, r, sigma, b)
    carry = math.exp((b - r) * t)
    disc = math.exp(-r * t)
    if ot is OptionType.CALL:
        return S * carry * norm_cdf(d1) - K * disc * norm_cdf(d2)
    return K * disc * norm_cdf(-d2) - S * carry * norm_cdf(-d1)


def call_price(S, K, t, r, sigma, b=None) -> float:
    return price(S, K, t, r, sigma, OptionType.CALL, b)


def put_price(S, K, t, r, sigma, b=None) -> float:
    return price(S, K, t, r, sigma, OptionType.PUT, b)


def delta(S, K, t, r, sigma, option_type=OptionType.CALL, b=None) -> float:
    """dPrice/dS."""
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if t == 0 or sigma == 0:
        fwd = S * math.exp(b * t)
        itm = fwd > K if ot is OptionType.CALL else fwd < K
        base = math.exp((b - r) * t) if itm else 0.0
        return base if ot is OptionType.CALL else -base
    d1, _ = _d1_d2(S, K, t, r, sigma, b)
    carry = math.exp((b - r) * t)
    if ot is OptionType.CALL:
        return carry * norm_cdf(d1)
    return carry * (norm_cdf(d1) - 1.0)


def gamma(S, K, t, r, sigma, b=None) -> float:
    """d2Price/dS2. Identical for calls and puts."""
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if t == 0 or sigma == 0:
        return 0.0
    d1, _ = _d1_d2(S, K, t, r, sigma, b)
    carry = math.exp((b - r) * t)
    return carry * norm_pdf(d1) / (S * sigma * math.sqrt(t))


def vega(S, K, t, r, sigma, b=None) -> float:
    """dPrice/dSigma, per 1.0 change in vol (divide by 100 for per-vol-point)."""
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if t == 0 or sigma == 0:
        return 0.0
    d1, _ = _d1_d2(S, K, t, r, sigma, b)
    carry = math.exp((b - r) * t)
    return S * carry * norm_pdf(d1) * math.sqrt(t)


def theta(S, K, t, r, sigma, option_type=OptionType.CALL, b=None) -> float:
    """Calendar-time theta, dPrice/d(calendar time) per year.

    This is the market convention: equal to ``-dPrice/dt_expiry``, so long
    options usually show negative theta (value decays as the clock advances).
    Divide by 365 for per-calendar-day decay.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if t == 0 or sigma == 0:
        return 0.0
    d1, d2 = _d1_d2(S, K, t, r, sigma, b)
    carry = math.exp((b - r) * t)
    disc = math.exp(-r * t)
    term1 = -S * carry * norm_pdf(d1) * sigma / (2.0 * math.sqrt(t))
    if ot is OptionType.CALL:
        term2 = -(b - r) * S * carry * norm_cdf(d1)
        term3 = -r * K * disc * norm_cdf(d2)
        return term1 + term2 + term3
    term2 = (b - r) * S * carry * norm_cdf(-d1)
    term3 = r * K * disc * norm_cdf(-d2)
    return term1 + term2 + term3


def rho(S, K, t, r, sigma, option_type=OptionType.CALL, b=None) -> float:
    """dPrice/dr, per 1.0 change in rate.

    Assumes carry moves with the rate (the plain BSM stock case). For models
    where ``b`` is fixed independently of ``r`` (e.g. Black-76), pass ``b`` and
    interpret accordingly.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if t == 0 or sigma == 0:
        return 0.0
    _, d2 = _d1_d2(S, K, t, r, sigma, b)
    disc = math.exp(-r * t)
    if ot is OptionType.CALL:
        return K * t * disc * norm_cdf(d2)
    return -K * t * disc * norm_cdf(-d2)


def epsilon(S, K, t, r, sigma, option_type=OptionType.CALL, b=None) -> float:
    """Dividend rho (a.k.a. epsilon / psi): dPrice/dq, per 1.0 change in the
    continuous dividend yield.

    The dividend yield enters through the carry ``b = r - q``, so raising q
    lowers the forward. For a call ``dPrice/dq = -S t e^{(b-r)t} N(d1)``; for a
    put ``+S t e^{(b-r)t} N(-d1)``. Assumes ``b`` moves with ``q`` (the standard
    dividend-yield case).
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if t == 0 or sigma == 0:
        return 0.0
    d1, _ = _d1_d2(S, K, t, r, sigma, b)
    carry = math.exp((b - r) * t)
    if ot is OptionType.CALL:
        return -S * t * carry * norm_cdf(d1)
    return S * t * carry * norm_cdf(-d1)


@dataclass(frozen=True)
class Greeks:
    price: float
    delta: float
    gamma: float
    vega: float
    theta: float
    rho: float


def greeks(S, K, t, r, sigma, option_type=OptionType.CALL, b=None) -> Greeks:
    """Compute price and all first/second-order Greeks in one call."""
    return Greeks(
        price=price(S, K, t, r, sigma, option_type, b),
        delta=delta(S, K, t, r, sigma, option_type, b),
        gamma=gamma(S, K, t, r, sigma, b),
        vega=vega(S, K, t, r, sigma, b),
        theta=theta(S, K, t, r, sigma, option_type, b),
        rho=rho(S, K, t, r, sigma, option_type, b),
    )
