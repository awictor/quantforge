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
