"""Chatterjee's xi rank correlation and Blomqvist's beta.

Two simple rank-based dependence measures:

  * ``chatterjee_xi`` -- Chatterjee's (2020) coefficient. It measures how much ``Y``
    is a *function* of ``X``: order the pairs by ``X``, look at how the ranks of the
    ``Y`` values jump between neighbours, and small jumps mean ``Y`` moves smoothly
    with ``X``. It is 0 under independence and tends to 1 when ``Y`` is a noiseless
    (possibly nonmonotone) function of ``X`` -- unlike Pearson/Kendall it is *not*
    symmetric and detects any functional relationship, oscillatory ones included.
  * ``blomqvist_beta`` -- the medial-correlation coefficient: the fraction of points
    in the two "concordant" quadrants around the medians, rescaled to ``[-1, 1]``.
    A fast, robust, median-based cousin of Kendall's tau.

Pure standard library.
"""

import math


def _median(values):
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2:
        return s[mid]
    return 0.5 * (s[mid - 1] + s[mid])


def chatterjee_xi(x, y):
    """Chatterjee's xi rank correlation of ``Y`` on ``X`` (asymmetric).

    Sorts by ``x`` (ties in ``x`` broken by stable order), then with
    ``r_i`` = #{j : y_j <= y_(i)} and ``l_i`` = #{j : y_j >= y_(i)},

        xi = 1 - n * sum_i |r_{i+1} - r_i| / (2 * sum_i l_i (n - l_i)).

    Zero under independence, approaching 1 when ``Y`` is a noiseless function of
    ``X``. Handles ties in ``y`` via the general (Azadkia-Chatterjee) form.
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have equal length")
    if n < 2:
        raise ValueError("need at least 2 points")

    order = sorted(range(n), key=lambda i: x[i])
    ys = [y[i] for i in order]

    # r_i = number of j with y_j <= ys[i]; l_i = number with y_j >= ys[i].
    r = [sum(1 for v in ys if v <= yi) for yi in ys]
    l = [sum(1 for v in ys if v >= yi) for yi in ys]

    num = n * sum(abs(r[i + 1] - r[i]) for i in range(n - 1))
    den = 2 * sum(li * (n - li) for li in l)
    if den == 0:
        return 0.0
    return 1.0 - num / den


def blomqvist_beta(x, y):
    """Blomqvist's beta (medial correlation).

    ``beta = (n_concordant - n_discordant) / n_used`` where a point is concordant if
    it sits in the same direction from both medians (both above or both below) and
    discordant otherwise; points exactly on a median are dropped (``n_used`` counts
    only the points kept). Lies in ``[-1, 1]``:
    +1 comonotone, -1 countermonotone, 0 under independence.
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have equal length")
    if n < 2:
        raise ValueError("need at least 2 points")
    mx = _median(x)
    my = _median(y)
    conc = disc = used = 0
    for i in range(n):
        sx = x[i] - mx
        sy = y[i] - my
        if sx == 0 or sy == 0:
            continue
        used += 1
        if (sx > 0) == (sy > 0):
            conc += 1
        else:
            disc += 1
    if used == 0:
        return 0.0
    return (conc - disc) / used
