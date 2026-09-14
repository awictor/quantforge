"""Latin hypercube sampling and uniformity diagnostics.

Latin hypercube sampling stratifies every dimension: split each axis into ``n`` equal
bins and place exactly one sample in each bin per axis, so the projection onto any single
coordinate is perfectly uniform. That removes the clustering and gaps of plain random
sampling with none of the structure of a grid, which is why it is the standard design for
Monte-Carlo experiments and sensitivity analysis. A ``centered`` variant puts each point
at its bin center; a ``maximin`` search picks the most-spread design over several tries.
The L2 star discrepancy measures how far a point set departs from perfect uniformity.
Pure standard library (deterministic LCG stream per seed).
"""

import math


def _lcg(seed):
    state = seed & 0x7FFFFFFF
    if state == 0:
        state = 1

    def _next():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return state / 0x80000000

    return _next


def latin_hypercube(n, dim, seed=1234567, centered=False):
    """Latin hypercube design of ``n`` points in ``[0, 1)^dim``.

    Each of the ``dim`` axes is split into ``n`` equal bins with exactly one sample per
    bin, then the per-axis bin assignments are independently permuted. With
    ``centered=True`` each sample sits at its bin center (deterministic given the
    permutation); otherwise it is placed uniformly at random within its bin. Returns a
    list of ``n`` points (each a length-``dim`` list).
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if dim < 1:
        raise ValueError("dim must be >= 1")
    rand = _lcg(seed)
    # For each dimension, a permutation of [0, n).
    cols = []
    for _ in range(dim):
        perm = list(range(n))
        # Fisher-Yates with the LCG stream.
        for i in range(n - 1, 0, -1):
            j = int(rand() * (i + 1))
            perm[i], perm[j] = perm[j], perm[i]
        cols.append(perm)
    points = []
    for i in range(n):
        pt = []
        for d in range(dim):
            bin_index = cols[d][i]
            offset = 0.5 if centered else rand()
            pt.append((bin_index + offset) / n)
        points.append(pt)
    return points


def maximin_lhs(n, dim, seed=1234567, tries=20):
    """Latin hypercube maximizing the minimum inter-point distance over ``tries`` designs.

    Generates ``tries`` Latin hypercubes (distinct seeds) and returns the one whose
    closest pair is farthest apart -- a more space-filling design than a single random
    LHS, at ``tries`` times the cost. Returns the chosen point list.
    """
    if tries < 1:
        raise ValueError("tries must be >= 1")
    best = None
    best_score = -1.0
    for t in range(tries):
        pts = latin_hypercube(n, dim, seed=seed + t * 100003)
        score = _min_pairwise_distance(pts)
        if score > best_score:
            best_score = score
            best = pts
    return best


def _min_pairwise_distance(points):
    n = len(points)
    if n < 2:
        return float("inf")
    best = float("inf")
    for i in range(n):
        for j in range(i + 1, n):
            d = math.sqrt(sum((points[i][k] - points[j][k]) ** 2
                              for k in range(len(points[i]))))
            if d < best:
                best = d
    return best


def l2_star_discrepancy(points):
    """L2 star discrepancy of a point set in ``[0, 1]^dim`` (Warnock's formula).

    A single number measuring departure from perfect uniformity -- smaller is more
    uniform. Uses Warnock's closed form, so it is exact and ``O(n^2 * dim)``. Useful for
    comparing a Latin hypercube (or Sobol/Halton set) against plain random sampling.
    """
    n = len(points)
    if n == 0:
        raise ValueError("need at least one point")
    d = len(points[0])
    term1 = 1.0 / (3.0 ** d)
    term2 = 0.0
    for p in points:
        prod = 1.0
        for x in p:
            prod *= (1.0 - x * x)
        term2 += prod
    term2 *= (2.0 ** (1 - d)) / n
    term3 = 0.0
    for pi in points:
        for pj in points:
            prod = 1.0
            for k in range(d):
                prod *= (1.0 - max(pi[k], pj[k]))
            term3 += prod
    term3 /= n * n
    return math.sqrt(abs(term1 - term2 + term3))
