"""Model-selection information criteria from a log-likelihood.

Generic, model-agnostic criteria that trade goodness-of-fit against complexity, so
the model minimizing one is preferred:

- ``aic``  -- Akaike: ``-2 logL + 2k``,
- ``aicc`` -- small-sample corrected AIC: ``AIC + 2k(k+1)/(n-k-1)``,
- ``bic``  -- Schwarz/Bayesian: ``-2 logL + k ln n`` (heavier penalty, favors
  parsimony),
- ``hqic`` -- Hannan-Quinn: ``-2 logL + 2k ln(ln n)``.

For a Gaussian model fit by least squares, ``gaussian_log_likelihood`` turns a
residual sum of squares into the log-likelihood these consume (the noise variance
counts as one estimated parameter). Pure standard library.
"""

import math


def gaussian_log_likelihood(rss, n):
    """Concentrated Gaussian log-likelihood from a residual sum of squares.

    ``logL = -n/2 (ln(2 pi) + ln(rss/n) + 1)`` at the MLE noise variance
    ``sigma^2 = rss/n``. Increases (toward zero from below) as the fit improves.
    """
    if n < 1:
        raise ValueError("n must be positive")
    if rss <= 0.0:
        raise ValueError("rss must be positive")
    sigma2 = rss / n
    return -0.5 * n * (math.log(2.0 * math.pi) + math.log(sigma2) + 1.0)


def aic(log_likelihood, k):
    """Akaike information criterion ``-2 logL + 2k`` (``k`` = number of parameters)."""
    if k < 0:
        raise ValueError("k must be non-negative")
    return -2.0 * log_likelihood + 2.0 * k


def aicc(log_likelihood, k, n):
    """Small-sample corrected AIC ``AIC + 2k(k+1)/(n-k-1)``.

    Approaches :func:`aic` as ``n`` grows; use it when ``n`` is not much larger than
    ``k``. Requires ``n > k + 1``.
    """
    if k < 0:
        raise ValueError("k must be non-negative")
    if n <= k + 1:
        raise ValueError("need n > k + 1 for the AICc correction")
    return aic(log_likelihood, k) + 2.0 * k * (k + 1) / (n - k - 1)


def bic(log_likelihood, k, n):
    """Bayesian (Schwarz) information criterion ``-2 logL + k ln n``.

    Penalizes complexity more than AIC for ``n >= 8`` (``ln n > 2``), so it selects
    more parsimonious models.
    """
    if k < 0:
        raise ValueError("k must be non-negative")
    if n < 1:
        raise ValueError("n must be positive")
    return -2.0 * log_likelihood + k * math.log(n)


def hqic(log_likelihood, k, n):
    """Hannan-Quinn information criterion ``-2 logL + 2k ln(ln n)``.

    A penalty between AIC and BIC for moderate ``n``. Requires ``n >= 3`` so that
    ``ln(ln n)`` is defined and positive.
    """
    if k < 0:
        raise ValueError("k must be non-negative")
    if n < 3:
        raise ValueError("need n >= 3 for the Hannan-Quinn penalty")
    return -2.0 * log_likelihood + 2.0 * k * math.log(math.log(n))
