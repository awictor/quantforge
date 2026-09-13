"""Two-sample Cramer-von Mises test of equal distributions.

The Cramer-von Mises criterion integrates the *squared* gap between two empirical
CDFs (where Kolmogorov-Smirnov takes only the single largest gap), so it uses the
whole curve and is more powerful against differences spread across the distribution.
Anderson's (1962) rank form of the two-sample statistic is

    T = U / (N M (N + M)) - (4 M N - 1) / (6 (N + M)),

where ``N``, ``M`` are the sample sizes, ``N + M`` the pool, and ``U`` a
rank-based sum. Under equal distributions ``T`` converges to the first
Cramer-von Mises limit distribution, whose upper tail gives the p-value. Univariate
samples; pure standard library.
"""

import math


def _ranks_in_pool(sample, pool_sorted):
    """Midranks of each value of ``sample`` within the combined sorted pool."""
    # Midrank handles ties: average of the 1-based positions of equal values.
    ranks = []
    n = len(pool_sorted)
    for v in sample:
        lo = 0
        hi = n
        # first index >= v
        while lo < hi:
            mid = (lo + hi) // 2
            if pool_sorted[mid] < v:
                lo = mid + 1
            else:
                hi = mid
        first = lo
        # last index <= v
        lo2, hi2 = 0, n
        while lo2 < hi2:
            mid = (lo2 + hi2) // 2
            if pool_sorted[mid] <= v:
                lo2 = mid + 1
            else:
                hi2 = mid
        last = lo2 - 1
        # 1-based midrank
        ranks.append((first + 1 + last + 1) / 2.0)
    return ranks


def _cvm_asymptotic_sf(t):
    """Upper tail P(T >= t) of the first limiting Cramer-von Mises distribution.

    Uses the standard series expansion (Anderson-Darling 1952) in terms of the
    modified Bessel function K_{1/4}, evaluated by a rapidly converging sum.
    """
    if t <= 0:
        return 1.0
    # Series from the CvM limit CDF; a handful of terms suffices in the usable range.
    # P(T < t) = (1/pi sqrt(t)) sum_{j>=0} (Gamma(j+1/2)/Gamma(1/2)/j!) sqrt(4j+1)
    #            exp(-(4j+1)^2/(16 t)) K_{1/4}((4j+1)^2/(16 t))
    total = 0.0
    for j in range(0, 100):
        c = math.gamma(j + 0.5) / (math.gamma(0.5) * math.factorial(j))
        a = (4 * j + 1) ** 2 / (16.0 * t)
        term = c * math.sqrt(4 * j + 1) * math.exp(-a) * _besselk_quarter(a)
        total += term
        if j > 2 and abs(term) < 1e-12 * (abs(total) + 1e-30):
            break
    cdf = total / (math.pi * math.sqrt(t))
    cdf = min(max(cdf, 0.0), 1.0)
    return 1.0 - cdf


def _besselk_quarter(x):
    """Modified Bessel function K_{1/4}(x) for x > 0 via numeric quadrature.

    Integral representation ``K_nu(x) = integral_0^inf exp(-x cosh t) cosh(nu t) dt``,
    evaluated by adaptive Simpson on a truncated range -- adequate for the moderate
    ``x`` arising in the CvM series.
    """
    nu = 0.25

    def f(u):
        return math.exp(-x * math.cosh(u)) * math.cosh(nu * u)

    # Integrand decays like exp(-x cosh u); truncate where it is negligible.
    upper = 1.0
    while f(upper) > 1e-16 and upper < 50.0:
        upper += 1.0
    n = 400
    h = upper / n
    s = f(0.0) + f(upper)
    for i in range(1, n):
        s += (4.0 if i % 2 else 2.0) * f(i * h)
    return s * h / 3.0


def cramer_von_mises_2samp(a, b):
    """Two-sample Cramer-von Mises test.

    Returns a dict with the ``statistic`` (Anderson's ``T``) and the asymptotic
    ``p_value`` from the limiting Cramer-von Mises distribution. A small p-value
    rejects the null that ``a`` and ``b`` are drawn from the same distribution.
    """
    n = len(a)
    m = len(b)
    if n == 0 or m == 0:
        raise ValueError("both samples must be non-empty")
    pool = sorted(list(a) + list(b))
    ra = _ranks_in_pool(a, pool)
    rb = _ranks_in_pool(b, pool)

    # Anderson (1962): U = N sum_i (r_i - i)^2 + M sum_j (s_j - j)^2, with the
    # sample values ordered and i, j their within-sample 1-based positions.
    a_sorted = sorted(range(n), key=lambda i: a[i])
    b_sorted = sorted(range(m), key=lambda j: b[j])
    ra_ord = [ra[i] for i in a_sorted]
    rb_ord = [rb[j] for j in b_sorted]
    u = n * sum((ra_ord[i] - (i + 1)) ** 2 for i in range(n)) \
        + m * sum((rb_ord[j] - (j + 1)) ** 2 for j in range(m))
    big_n = n + m
    t = u / (n * m * big_n) - (4 * m * n - 1) / (6.0 * big_n)
    p = _cvm_asymptotic_sf(t)
    return {"statistic": t, "p_value": p}
