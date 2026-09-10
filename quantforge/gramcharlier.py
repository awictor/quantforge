"""Corrado-Su skewness/kurtosis-adjusted option pricing (Gram-Charlier).

Black-Scholes assumes lognormal returns (zero skew, zero excess kurtosis). Real
return distributions are skewed and fat-tailed. Corrado and Su (1996) add the
first skewness and kurtosis corrections to the Black-Scholes call via a
Gram-Charlier Type A expansion of the terminal density:

    C = C_BS + skew * Q3 + (kurt_excess) * Q4,

where ``Q3`` and ``Q4`` are closed-form adjustment terms built from the normal
pdf. Setting skew = 0 and excess kurtosis = 0 recovers Black-Scholes exactly.
Puts follow from put-call parity. Also included are the realized skewness and
excess-kurtosis estimators for a return series, so a historical distribution
can feed the price directly.
"""

import math
from typing import Sequence

from .mathfns import norm_cdf, norm_pdf
from .bsm import call_price, OptionType, _coerce_type, _validate


def corrado_su_call(S, K, t, r, sigma, skew=0.0, excess_kurt=0.0, b=None) -> float:
    """Corrado-Su (1996) skew/kurtosis-adjusted European call price.

    Args:
        skew: skewness of the (log) return distribution.
        excess_kurt: excess kurtosis (kurtosis - 3).
        b: cost of carry (defaults to r). skew=kurt=0 => Black-Scholes.
    """
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    bs = call_price(S, K, t, r, sigma, b=b)
    if t == 0 or sigma == 0:
        return bs

    vsqrt = sigma * math.sqrt(t)
    d = (math.log(S / K) + (b + 0.5 * sigma * sigma) * t) / vsqrt
    carry = math.exp((b - r) * t)
    nd = norm_pdf(d)

    # Corrado-Su adjustment terms (using the "corrected" 1997 coefficients).
    Q3 = (1.0 / math.factorial(3)) * S * carry * vsqrt * (
        (2.0 * vsqrt - d) * nd)
    Q4 = (1.0 / math.factorial(4)) * S * carry * vsqrt * (
        (d * d - 1.0 - 3.0 * vsqrt * (d - vsqrt)) * nd)

    return bs + skew * Q3 + excess_kurt * Q4


def corrado_su_price(S, K, t, r, sigma, skew=0.0, excess_kurt=0.0,
                     option_type=OptionType.CALL, b=None) -> float:
    """Corrado-Su price for a call or put (put via put-call parity)."""
    ot = _coerce_type(option_type)
    if b is None:
        b = r
    call = corrado_su_call(S, K, t, r, sigma, skew, excess_kurt, b)
    if ot is OptionType.CALL:
        return call
    return call - S * math.exp((b - r) * t) + K * math.exp(-r * t)


def realized_skewness(returns: Sequence[float]) -> float:
    """Sample skewness of a return series (bias-corrected denominator n)."""
    n = len(returns)
    if n < 3:
        raise ValueError("need at least three returns")
    mean = sum(returns) / n
    m2 = sum((x - mean) ** 2 for x in returns) / n
    m3 = sum((x - mean) ** 3 for x in returns) / n
    if m2 <= 0:
        return 0.0
    return m3 / (m2 ** 1.5)


def realized_excess_kurtosis(returns: Sequence[float]) -> float:
    """Sample excess kurtosis (kurtosis - 3) of a return series."""
    n = len(returns)
    if n < 4:
        raise ValueError("need at least four returns")
    mean = sum(returns) / n
    m2 = sum((x - mean) ** 2 for x in returns) / n
    m4 = sum((x - mean) ** 4 for x in returns) / n
    if m2 <= 0:
        return 0.0
    return m4 / (m2 * m2) - 3.0
