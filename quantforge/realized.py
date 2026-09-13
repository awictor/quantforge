"""Jump-robust realized volatility measures from intraday returns.

Given a series of high-frequency (log) returns within a period, the realized
variance ``RV = sum r_i^2`` consistently estimates the total quadratic variation
-- but it lumps the continuous diffusion together with any price jumps. Barndorff-
Nielsen and Shephard's bipower variation

    BV = (pi / 2) * sum_{i=2}^{n} |r_i| |r_{i-1}|

converges instead to the *integrated variance* alone and is robust to a finite
number of jumps (a jump enters only one product term and is washed out
asymptotically). Their difference isolates the jump contribution:

    jump variation = max(RV - BV, 0).

These are the standard building blocks of realized-volatility jump tests. Pure
standard library.
"""

import math

_MU1 = math.sqrt(2.0 / math.pi)         # E|Z| for standard normal Z
_BV_SCALE = 1.0 / (_MU1 * _MU1)         # = pi / 2


def realized_variance_from_returns(returns):
    """Realized variance ``sum r_i^2`` of a return series.

    Consistent for the total quadratic variation (diffusion plus jumps) as the
    sampling frequency rises.
    """
    if len(returns) == 0:
        raise ValueError("need at least 1 return")
    return sum(r * r for r in returns)


def bipower_variation(returns):
    """Barndorff-Nielsen-Shephard bipower variation.

    ``BV = (pi/2) sum_{i=2}^{n} |r_i| |r_{i-1}|``. Estimates the integrated
    variance of the continuous part only and is robust to jumps. Requires at
    least two returns.
    """
    n = len(returns)
    if n < 2:
        raise ValueError("need at least 2 returns")
    s = 0.0
    for i in range(1, n):
        s += abs(returns[i]) * abs(returns[i - 1])
    return _BV_SCALE * s


def jump_variation(returns):
    """Jump component of the quadratic variation, ``max(RV - BV, 0)``.

    Zero (up to sampling noise) for a purely continuous path; strictly positive
    when the returns contain a jump. Clamped at zero because the estimator can go
    slightly negative on jump-free data.
    """
    return max(realized_variance_from_returns(returns) - bipower_variation(returns),
               0.0)


def realized_volatility_signature(returns, annualization=1.0):
    """Realized volatility ``sqrt(RV)``, optionally annualized.

    ``annualization`` multiplies the variance before the square root (e.g. the
    number of periods per year for intraday returns aggregated to one day times
    252). Defaults to 1 (the raw realized vol of the supplied returns).
    """
    return math.sqrt(realized_variance_from_returns(returns) * annualization)


# Andersen-Dobrev-Schaumburg (2012) nearest-neighbour truncation constants.
_MINRV_SCALE = math.pi / (math.pi - 2.0)
_MEDRV_SCALE = math.pi / (6.0 - 4.0 * math.sqrt(3.0) + math.pi)


def min_realized_variance(returns):
    """MinRV jump-robust integrated-variance estimator (Andersen-Dobrev-Schaumburg).

    ``MinRV = (pi / (pi - 2)) * (n / (n - 1)) * sum_i min(|r_i|, |r_{i+1}|)^2``. Each
    term pairs adjacent returns and keeps the smaller magnitude, so an isolated jump
    (which lands in one return) is discarded by the minimum. Converges to the
    integrated variance of the continuous part; more robust to jumps than bipower
    variation and to occasional zero returns. Requires at least two returns.
    """
    n = len(returns)
    if n < 2:
        raise ValueError("need at least 2 returns")
    s = 0.0
    for i in range(n - 1):
        a, b = abs(returns[i]), abs(returns[i + 1])
        m = a if a < b else b
        s += m * m
    return _MINRV_SCALE * (n / (n - 1.0)) * s


def med_realized_variance(returns):
    """MedRV jump-robust integrated-variance estimator (Andersen-Dobrev-Schaumburg).

    ``MedRV = c * (n / (n - 2)) * sum_i median(|r_{i-1}|, |r_i|, |r_{i+1}|)^2`` with
    ``c = pi / (6 - 4 sqrt(3) + pi)``. Taking the median of three neighbouring
    magnitudes discards a lone jump and, unlike MinRV, is also robust to two nearby
    jumps and less sensitive to zero returns. Requires at least three returns.
    """
    n = len(returns)
    if n < 3:
        raise ValueError("need at least 3 returns")
    s = 0.0
    for i in range(1, n - 1):
        trio = sorted((abs(returns[i - 1]), abs(returns[i]), abs(returns[i + 1])))
        med = trio[1]
        s += med * med
    return _MEDRV_SCALE * (n / (n - 2.0)) * s


def realized_quarticity(returns):
    """Realized quarticity ``(n / 3) * sum r_i^4``.

    A consistent estimator of the integrated quarticity ``integral sigma^4``, which
    sets the asymptotic variance of realized variance and appears in the standard
    errors of realized-volatility jump tests. Non-negative.
    """
    n = len(returns)
    if n == 0:
        raise ValueError("need at least 1 return")
    return (n / 3.0) * sum(r ** 4 for r in returns)
