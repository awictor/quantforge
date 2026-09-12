"""Chain-ladder loss reserving from a claims development triangle.

The chain-ladder method projects a cumulative-claims run-off triangle to ultimate
losses using volume-weighted age-to-age development factors:

    f_j = sum_i C_{i, j+1} / sum_i C_{i, j}   (over rows with both entries),

then fills each row forward, ``C_{i, k+1} = C_{i, k} * f_k``. The reserve (IBNR)
for an accident year is ``ultimate - latest paid``. Standard non-life reserving.
The triangle is a list of rows; row ``i`` has ``n - i`` observed cumulative
entries (most recent accident year is the shortest). Pure standard library.
"""


def development_factors(triangle):
    """Volume-weighted age-to-age development factors from a cumulative triangle.

    ``triangle[i]`` is the observed cumulative claims for accident year ``i`` at
    development ages ``0 .. len(triangle[i]) - 1``. Returns ``n - 1`` factors
    ``f_0 .. f_{n-2}`` linking successive development ages.
    """
    n = len(triangle)
    if n < 2:
        raise ValueError("need at least 2 accident years")
    factors = []
    for j in range(n - 1):
        num = 0.0
        den = 0.0
        for i in range(n):
            if len(triangle[i]) > j + 1:
                den += triangle[i][j]
                num += triangle[i][j + 1]
        if den <= 0.0:
            raise ValueError("no data to estimate a development factor")
        factors.append(num / den)
    return factors


def chain_ladder(triangle):
    """Project a claims triangle to ultimate losses and reserves.

    Returns a dict with ``factors`` (age-to-age), ``ultimate`` (per accident year),
    ``reserve`` (IBNR per year = ultimate - latest observed), and
    ``total_reserve``. A fully-developed row has zero reserve.
    """
    n = len(triangle)
    f = development_factors(triangle)
    ultimate = []
    reserve = []
    for i in range(n):
        row = list(triangle[i])
        latest = row[-1]
        age = len(row) - 1
        proj = latest
        for k in range(age, n - 1):
            proj *= f[k]
        ultimate.append(proj)
        reserve.append(proj - latest)
    return {
        "factors": f,
        "ultimate": ultimate,
        "reserve": reserve,
        "total_reserve": sum(reserve),
    }


def development_pattern(factors):
    """Cumulative development pattern (% reported) from age-to-age factors.

    Returns ``pct[j]`` = fraction of ultimate developed by age ``j``, computed as
    the reciprocal of the cumulative product of the remaining factors. The final
    age is fully developed (``1.0``).
    """
    n = len(factors) + 1
    pct = [0.0] * n
    pct[n - 1] = 1.0
    cdf = 1.0
    for j in range(n - 2, -1, -1):
        cdf *= factors[j]
        pct[j] = 1.0 / cdf
    return pct


def bornhuetter_ferguson(triangle, apriori_ultimates):
    """Bornhuetter-Ferguson reserves blending development with a-priori ultimates.

    For each accident year the BF reserve is
    ``apriori_ultimate * (1 - pct_developed)``, where ``pct_developed`` comes from
    the chain-ladder development pattern. The BF ultimate is the latest paid plus
    that reserve -- an a-priori-anchored estimate that is robust for green
    (little-developed) years where chain-ladder is volatile.

    Parameters
    ----------
    triangle : list[list[float]]
        Cumulative-claims triangle (as in :func:`chain_ladder`).
    apriori_ultimates : sequence of float
        A-priori ultimate loss per accident year (e.g. premium x expected loss
        ratio).

    Returns
    -------
    dict
        ``pattern`` (% developed by age), ``reserve`` and ``ultimate`` per year,
        and ``total_reserve``.
    """
    n = len(triangle)
    if len(apriori_ultimates) != n:
        raise ValueError("apriori_ultimates must match the number of accident years")
    f = development_factors(triangle)
    pattern = development_pattern(f)
    reserve = []
    ultimate = []
    for i in range(n):
        age = len(triangle[i]) - 1
        pct = pattern[age]
        res = apriori_ultimates[i] * (1.0 - pct)
        reserve.append(res)
        ultimate.append(triangle[i][-1] + res)
    return {
        "pattern": pattern,
        "reserve": reserve,
        "ultimate": ultimate,
        "total_reserve": sum(reserve),
    }


def cape_cod(triangle, premiums):
    """Cape Cod (Stanard-Buhlmann) reserving.

    Like Bornhuetter-Ferguson but the a-priori loss ratio is estimated from the
    data rather than assumed: the expected loss ratio is

        ELR = sum_i latest_i / sum_i (premium_i * pct_developed_i),

    the total observed losses over the total "used-up" premium (premium weighted
    by how developed each year is). Each year's a-priori ultimate is then
    ``premium_i * ELR`` and its reserve ``apriori * (1 - pct_developed_i)``.

    Parameters
    ----------
    triangle : list[list[float]]
        Cumulative-claims triangle (as in :func:`chain_ladder`).
    premiums : sequence of float
        Earned premium per accident year.

    Returns
    -------
    dict
        ``elr``, ``pattern``, ``reserve`` / ``ultimate`` per year, and
        ``total_reserve``.
    """
    n = len(triangle)
    if len(premiums) != n:
        raise ValueError("premiums must match the number of accident years")
    f = development_factors(triangle)
    pattern = development_pattern(f)
    pct = [pattern[len(triangle[i]) - 1] for i in range(n)]

    total_loss = sum(triangle[i][-1] for i in range(n))
    used_premium = sum(premiums[i] * pct[i] for i in range(n))
    if used_premium <= 0.0:
        raise ValueError("total used-up premium must be positive")
    elr = total_loss / used_premium

    reserve = []
    ultimate = []
    for i in range(n):
        apriori = premiums[i] * elr
        res = apriori * (1.0 - pct[i])
        reserve.append(res)
        ultimate.append(triangle[i][-1] + res)
    return {
        "elr": elr,
        "pattern": pattern,
        "reserve": reserve,
        "ultimate": ultimate,
        "total_reserve": sum(reserve),
    }


def incremental_to_cumulative(triangle):
    """Convert an incremental-claims triangle to cumulative.

    Each row's cumulative entry is the running sum of its incremental entries.
    Ragged rows (shorter for recent accident years) are preserved.
    """
    out = []
    for row in triangle:
        cum = []
        running = 0.0
        for v in row:
            running += v
            cum.append(running)
        out.append(cum)
    return out


def cumulative_to_incremental(triangle):
    """Convert a cumulative-claims triangle to incremental (successive differences)."""
    out = []
    for row in triangle:
        inc = []
        prev = 0.0
        for v in row:
            inc.append(v - prev)
            prev = v
        out.append(inc)
    return out


def paid_to_date(cumulative_triangle):
    """Latest (diagonal) paid amount per accident year of a cumulative triangle."""
    return [row[-1] for row in cumulative_triangle]
