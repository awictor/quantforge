"""One-dimensional Wasserstein (earth-mover) distance between samples.

The p-Wasserstein distance is the minimum "work" to morph one distribution into
another, moving probability mass over the ground distance. In 1-D it has a closed form
via the *quantile functions*: for equal-size samples the p-Wasserstein distance is just
the p-norm of the gap between the two sorted samples,

    W_p(x, y) = ( mean_i |x_(i) - y_(i)|^p )^{1/p}.

For unequal sizes (or weighted distributions) it integrates ``|F_x^{-1} - F_y^{-1}|^p``
over ``[0, 1]`` by merging the two sample CDFs. Unlike KL it is a true metric, finite
even for disjoint supports, and reflects *how far* mass moved -- ideal for comparing
return distributions or histograms. Pure standard library.
"""


def wasserstein1_sorted(x, y):
    """1-Wasserstein distance for two equal-length samples (sorted-difference form).

    ``W_1 = mean_i |x_(i) - y_(i)|`` after sorting each sample. Requires equal lengths;
    use :func:`wasserstein_distance` for unequal sizes.
    """
    if len(x) != len(y):
        raise ValueError("samples must have equal length; use wasserstein_distance")
    if not x:
        raise ValueError("need at least one point")
    xs = sorted(x)
    ys = sorted(y)
    return sum(abs(xs[i] - ys[i]) for i in range(len(xs))) / len(xs)


def wasserstein_distance(x, y, p=1):
    """p-Wasserstein distance between two 1-D samples of any sizes.

    Merges the empirical CDFs and integrates ``|F_x^{-1}(u) - F_y^{-1}(u)|^p`` over
    ``u in [0, 1]`` via the standard interval decomposition (SciPy's approach). ``p=1``
    is earth-mover distance; ``p=2`` the quadratic transport cost. Returns a true metric
    for ``p=1``. Both samples must be non-empty.
    """
    if not x or not y:
        raise ValueError("both samples must be non-empty")
    if p < 1:
        raise ValueError("p must be >= 1")
    xs = sorted(x)
    ys = sorted(y)
    nx, ny = len(xs), len(ys)

    # Merge the sorted CDF breakpoints; between consecutive u-levels each quantile is
    # constant, so accumulate |dx - dy|^p * (width of u-interval).
    all_vals = sorted(set(xs) | set(ys))
    # Cumulative-probability position of each sample value.
    total = 0.0
    ix = iy = 0
    prev_cdf = 0.0
    # Walk through combined support; between adjacent points the u-interval carries
    # the current quantile difference. Use the CDF-step decomposition.
    deltas = _merged_deltas(xs, ys)
    for width, xq, yq in deltas:
        total += (abs(xq - yq) ** p) * width
    return total ** (1.0 / p)


def _merged_deltas(xs, ys):
    """Yield (u-width, x-quantile, y-quantile) over the merged CDF grid."""
    nx, ny = len(xs), len(ys)
    # All CDF jump levels from both empirical distributions, in [0, 1].
    levels = sorted(set([i / nx for i in range(nx + 1)] + [j / ny for j in range(ny + 1)]))
    out = []
    for k in range(len(levels) - 1):
        lo, hi = levels[k], levels[k + 1]
        mid = 0.5 * (lo + hi)
        # Quantile at u=mid: value at index floor(mid*n) (clamped).
        xi = min(int(mid * nx), nx - 1)
        yi = min(int(mid * ny), ny - 1)
        out.append((hi - lo, xs[xi], ys[yi]))
    return out
