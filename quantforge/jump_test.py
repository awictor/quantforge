"""Barndorff-Nielsen-Shephard / Huang-Tauchen realized-volatility jump test.

Tests the null of *no jump* in a day of high-frequency returns. Realized variance
``RV`` captures the total quadratic variation while bipower variation ``BV`` captures
only the continuous part, so the relative jump ``(RV - BV) / RV`` is centred at zero
under the null. Scaled by its asymptotic standard error -- built from the tripower
quarticity ``TQ`` -- the ratio statistic

    z = sqrt(n) * (RV - BV) / RV / sqrt((mu1^{-4} + 2 mu1^{-2} - 5) * max(1, TQ / BV^2))

is asymptotically standard normal. A large positive ``z`` (small upper-tail
p-value) rejects "no jump". Pure standard library on top of :mod:`quantforge.realized`.
"""

import math

from .realized import realized_variance_from_returns, bipower_variation
from .mathfns import norm_cdf

_MU1 = math.sqrt(2.0 / math.pi)          # E|Z|
_THETA = math.pi * math.pi / 4.0 + math.pi - 5.0    # = mu1^{-4} + 2 mu1^{-2} - 5


def tripower_quarticity(returns):
    """Tripower quarticity ``TQ``, a jump-robust estimator of integrated quarticity.

    ``TQ = n * mu_{4/3}^{-3} * sum |r_{i-2}|^{4/3} |r_{i-1}|^{4/3} |r_i|^{4/3}`` with
    ``mu_{4/3} = 2^{2/3} Gamma(7/6) / Gamma(1/2)``. Robust to jumps (each enters only
    one triple), it estimates ``integral sigma^4`` and sets the scale of the jump
    test. Requires at least three returns.
    """
    n = len(returns)
    if n < 3:
        raise ValueError("need at least 3 returns")
    mu = 2.0 ** (2.0 / 3.0) * math.gamma(7.0 / 6.0) / math.gamma(0.5)
    p = 4.0 / 3.0
    s = 0.0
    for i in range(2, n):
        s += (abs(returns[i - 2]) ** p * abs(returns[i - 1]) ** p
              * abs(returns[i]) ** p)
    return n * (mu ** -3) * s


def bns_jump_test(returns):
    """Barndorff-Nielsen-Shephard ratio jump-test statistic and p-value.

    Returns ``(z, p_value)`` where ``z`` is the ratio statistic (asymptotically
    ``N(0, 1)`` under no jump) and ``p_value`` is the upper-tail ``P(Z > z)``. A
    small p-value rejects the no-jump null in favour of a jump having occurred
    during the sampled period. Requires at least three returns.
    """
    n = len(returns)
    if n < 3:
        raise ValueError("need at least 3 returns")
    rv = realized_variance_from_returns(returns)
    bv = bipower_variation(returns)
    if rv <= 0.0:
        raise ValueError("realized variance must be positive")
    tq = tripower_quarticity(returns)
    rel_jump = (rv - bv) / rv
    denom = _THETA * max(1.0, tq / (bv * bv)) / n if bv > 0.0 else _THETA / n
    z = rel_jump / math.sqrt(denom)
    return z, 1.0 - norm_cdf(z)
