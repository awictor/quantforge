"""Student's t distribution and fat-tailed parametric VaR / expected shortfall.

The Student-t captures the heavy tails of financial returns better than the
normal. This module provides the t density, CDF, and quantile (inverse CDF), and
the parametric value-at-risk and expected shortfall of a location-scale t. As the
degrees of freedom grow the t converges to the normal. Pure standard library.
"""

import math

from .mathfns import norm_cdf, norm_ppf


def _lgamma(x):
    return math.lgamma(x)


def t_pdf(x, df):
    """Student-t probability density with ``df`` degrees of freedom.

    ``f(x) = Gamma((df+1)/2) / (sqrt(df pi) Gamma(df/2)) (1 + x^2/df)^{-(df+1)/2}``.
    Symmetric about zero, heavier-tailed than the normal for finite ``df``.
    """
    if df <= 0:
        raise ValueError("df must be positive")
    c = math.exp(_lgamma((df + 1) / 2.0) - _lgamma(df / 2.0)) \
        / math.sqrt(df * math.pi)
    return c * (1.0 + x * x / df) ** (-(df + 1) / 2.0)


def t_cdf(x, df):
    """Student-t cumulative distribution via the regularized incomplete beta.

    Uses the identity ``P(T <= x) = 1 - 0.5 I_{df/(df+x^2)}(df/2, 1/2)`` for
    ``x > 0`` and symmetry for ``x < 0``. Converges to :func:`norm_cdf` as
    ``df -> inf``.
    """
    if df <= 0:
        raise ValueError("df must be positive")
    if x == 0.0:
        return 0.5
    xt = df / (df + x * x)
    ib = 0.5 * _betainc(df / 2.0, 0.5, xt)
    return 1.0 - ib if x > 0 else ib


def _betainc(a, b, x):
    """Regularized incomplete beta I_x(a, b) via the continued fraction."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = _lgamma(a) + _lgamma(b) - _lgamma(a + b)
    front = math.exp(math.log(x) * a + math.log(1.0 - x) * b - lbeta) / a
    # Lentz's continued fraction.
    f, c, d = 1.0, 1.0, 0.0
    for i in range(0, 200):
        m = i // 2
        if i == 0:
            num = 1.0
        elif i % 2 == 0:
            num = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
        else:
            num = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1.0 + num * d
        if abs(d) < 1e-30:
            d = 1e-30
        d = 1.0 / d
        c = 1.0 + num / c
        if abs(c) < 1e-30:
            c = 1e-30
        f *= c * d
        if abs(1.0 - c * d) < 1e-12:
            break
    result = front * (f - 1.0)
    # I_x(a,b); use symmetry if converging poorly for x > (a+1)/(a+b+2).
    if x < (a + 1.0) / (a + b + 2.0):
        return result
    return 1.0 - _betainc(b, a, 1.0 - x)


def t_ppf(p, df):
    """Student-t quantile (inverse CDF) by bisection on :func:`t_cdf`.

    Returns the ``x`` with ``t_cdf(x, df) = p``. Symmetric: ``t_ppf(1-p) =
    -t_ppf(p)``.
    """
    if not (0.0 < p < 1.0):
        raise ValueError("p must be in (0, 1)")
    if df <= 0:
        raise ValueError("df must be positive")
    lo, hi = -1e4, 1e4
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if t_cdf(mid, df) < p:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-10:
            break
    return 0.5 * (lo + hi)


def fit_df_from_kurtosis(excess_kurtosis):
    """Degrees of freedom implied by a sample's excess kurtosis (moment match).

    For a Student-t the excess kurtosis is ``6 / (df - 4)`` (finite only for
    ``df > 4``), so ``df = 4 + 6 / excess_kurtosis``. Requires positive excess
    kurtosis (fatter than normal); larger kurtosis implies fewer degrees of
    freedom (heavier tails).
    """
    if excess_kurtosis <= 0.0:
        raise ValueError("excess_kurtosis must be positive (fatter than normal)")
    return 4.0 + 6.0 / excess_kurtosis


def fit_student_t(returns):
    """Fit a location-scale Student-t to a return sample by moment matching.

    Matches the sample mean, variance, and excess kurtosis: ``df`` from
    :func:`fit_df_from_kurtosis`, and the scale from ``variance = scale^2 df/(df-2)``.
    Returns ``(mean, scale, df)``. Requires ``df > 4`` (positive excess kurtosis).
    """
    n = len(returns)
    if n < 4:
        raise ValueError("need at least four observations")
    mean = sum(returns) / n
    var = sum((x - mean) ** 2 for x in returns) / n
    if var <= 0.0:
        raise ValueError("zero-variance returns")
    m4 = sum((x - mean) ** 4 for x in returns) / n
    excess = m4 / (var * var) - 3.0
    df = fit_df_from_kurtosis(excess)
    if df <= 2.0:
        raise ValueError("implied df <= 2; variance is undefined for the t")
    scale = math.sqrt(var * (df - 2.0) / df)
    return mean, scale, df


def student_t_var(mean, scale, df, confidence=0.95):
    """Parametric VaR of a location-scale Student-t (positive loss magnitude).

    ``-(mean + scale * t_ppf(1 - confidence, df))`` -- heavier-tailed than the
    normal VaR for finite ``df``, converging to it as ``df -> inf``.
    """
    if scale <= 0:
        raise ValueError("scale must be positive")
    q = t_ppf(1.0 - confidence, df)
    return -(mean + scale * q)


def student_t_expected_shortfall(mean, scale, df, confidence=0.95):
    """Parametric expected shortfall of a location-scale Student-t.

    Closed form ``ES = -mean + scale * (df + q^2)/(df - 1) * f(q)/(1 - confidence)``
    where ``q = t_ppf(1 - confidence, df)`` and ``f`` is the t density. Requires
    ``df > 1`` (finite mean); always at least the :func:`student_t_var`.
    """
    if scale <= 0:
        raise ValueError("scale must be positive")
    if df <= 1:
        raise ValueError("df must exceed 1 for a finite expected shortfall")
    alpha = 1.0 - confidence
    q = t_ppf(alpha, df)
    tail = (df + q * q) / (df - 1.0) * t_pdf(q, df) / alpha
    return -mean + scale * tail
