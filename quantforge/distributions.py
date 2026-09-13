"""Common probability distributions built on the special functions.

CDFs, PDFs/PMFs and quantiles for the gamma, chi-square, Poisson, F and binomial
distributions, expressed through the regularized incomplete gamma and beta
functions in :mod:`quantforge.special`. Each continuous CDF is monotone in its
argument, so the quantile is found by bisection on the CDF. Pure standard library.
"""

import math

from .special import gammainc, gammaincc, betainc


def gamma_cdf(x, shape, scale=1.0):
    """CDF of the gamma distribution ``Gamma(shape, scale)`` at ``x``.

    ``F(x) = P(shape, x / scale)`` via the regularized lower incomplete gamma.
    Reduces to the exponential CDF ``1 - e^{-x/scale}`` when ``shape = 1``.
    """
    if shape <= 0.0 or scale <= 0.0:
        raise ValueError("shape and scale must be positive")
    if x <= 0.0:
        return 0.0
    return gammainc(shape, x / scale)


def gamma_pdf(x, shape, scale=1.0):
    """Density of the gamma distribution at ``x >= 0``."""
    if shape <= 0.0 or scale <= 0.0:
        raise ValueError("shape and scale must be positive")
    if x < 0.0:
        return 0.0
    if x == 0.0:
        return 0.0 if shape > 1.0 else (1.0 / scale if shape == 1.0 else math.inf)
    return math.exp((shape - 1.0) * math.log(x) - x / scale
                    - shape * math.log(scale) - math.lgamma(shape))


def _bisect_quantile(cdf, p, lo, hi):
    """Bisection inverse of a monotone CDF, expanding ``hi`` until it brackets ``p``."""
    while cdf(hi) < p:
        hi *= 2.0
        if hi > 1e300:
            return hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if cdf(mid) < p:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-12 * (1.0 + abs(mid)):
            break
    return 0.5 * (lo + hi)


def gamma_ppf(p, shape, scale=1.0):
    """Quantile (inverse CDF) of the gamma distribution for ``p`` in ``(0, 1)``."""
    if not (0.0 < p < 1.0):
        raise ValueError("p must be in (0, 1)")
    if shape <= 0.0 or scale <= 0.0:
        raise ValueError("shape and scale must be positive")
    return _bisect_quantile(lambda x: gamma_cdf(x, shape, scale), p, 0.0,
                            max(1.0, shape) * scale)


def chi2_cdf(x, df):
    """CDF of the chi-square distribution with ``df`` degrees of freedom.

    A gamma with ``shape = df/2`` and ``scale = 2``: ``F(x) = P(df/2, x/2)``.
    """
    if df <= 0.0:
        raise ValueError("df must be positive")
    return gamma_cdf(x, df / 2.0, 2.0)


def chi2_sf(x, df):
    """Survival function ``1 - chi2_cdf(x, df)`` (upper tail, for p-values)."""
    if df <= 0.0:
        raise ValueError("df must be positive")
    if x <= 0.0:
        return 1.0
    return gammaincc(df / 2.0, x / 2.0)


def chi2_ppf(p, df):
    """Quantile of the chi-square distribution for ``p`` in ``(0, 1)``."""
    if not (0.0 < p < 1.0):
        raise ValueError("p must be in (0, 1)")
    if df <= 0.0:
        raise ValueError("df must be positive")
    return gamma_ppf(p, df / 2.0, 2.0)


def poisson_pmf(k, lam):
    """Poisson probability mass ``P(N = k) = e^{-lam} lam^k / k!``."""
    if lam < 0.0:
        raise ValueError("lam must be non-negative")
    if k < 0:
        return 0.0
    return math.exp(k * math.log(lam) - lam - math.lgamma(k + 1)) if lam > 0.0 else (
        1.0 if k == 0 else 0.0)


def poisson_cdf(k, lam):
    """Poisson CDF ``P(N <= k)`` via the gamma relation ``= Q(k+1, lam)``.

    Uses ``P(N <= k) = gammaincc(k + 1, lam)`` (the regularized upper incomplete
    gamma), exact and stable for large ``lam`` where summing masses would lose
    precision. ``k`` is floored to an integer.
    """
    if lam < 0.0:
        raise ValueError("lam must be non-negative")
    kk = math.floor(k)
    if kk < 0:
        return 0.0
    if lam == 0.0:
        return 1.0
    return gammaincc(kk + 1.0, lam)


def f_cdf(x, d1, d2):
    """CDF of the F distribution with ``(d1, d2)`` degrees of freedom.

    ``F(x) = I_{d1 x / (d1 x + d2)}(d1/2, d2/2)`` via the regularized incomplete
    beta.
    """
    if d1 <= 0.0 or d2 <= 0.0:
        raise ValueError("degrees of freedom must be positive")
    if x <= 0.0:
        return 0.0
    y = d1 * x / (d1 * x + d2)
    return betainc(d1 / 2.0, d2 / 2.0, y)


def f_ppf(p, d1, d2):
    """Quantile of the F distribution for ``p`` in ``(0, 1)``."""
    if not (0.0 < p < 1.0):
        raise ValueError("p must be in (0, 1)")
    if d1 <= 0.0 or d2 <= 0.0:
        raise ValueError("degrees of freedom must be positive")
    return _bisect_quantile(lambda x: f_cdf(x, d1, d2), p, 0.0, 2.0)


def binomial_cdf(k, n, prob):
    """Binomial CDF ``P(X <= k)`` for ``n`` trials with success probability ``prob``.

    Uses the incomplete-beta identity ``P(X <= k) = I_{1-p}(n - k, k + 1)``, exact
    and stable for large ``n``. ``k`` is floored to an integer.
    """
    if n < 0:
        raise ValueError("n must be non-negative")
    if not (0.0 <= prob <= 1.0):
        raise ValueError("prob must be in [0, 1]")
    kk = math.floor(k)
    if kk < 0:
        return 0.0
    if kk >= n:
        return 1.0
    return betainc(n - kk, kk + 1.0, 1.0 - prob)


def binomial_pmf(k, n, prob):
    """Binomial probability mass ``C(n, k) p^k (1-p)^{n-k}``."""
    if n < 0:
        raise ValueError("n must be non-negative")
    if not (0.0 <= prob <= 1.0):
        raise ValueError("prob must be in [0, 1]")
    if k < 0 or k > n:
        return 0.0
    log_coef = math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    if prob == 0.0:
        return 1.0 if k == 0 else 0.0
    if prob == 1.0:
        return 1.0 if k == n else 0.0
    return math.exp(log_coef + k * math.log(prob) + (n - k) * math.log(1.0 - prob))
