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
