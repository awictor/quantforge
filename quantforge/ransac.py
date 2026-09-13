"""RANSAC robust line fitting (random sample consensus).

When outliers exceed even the 50% breakdown of median estimators, RANSAC still finds
the true line: repeatedly fit a candidate from a *minimal* random sample (two points),
count how many of all points fall within a residual threshold of it (the *inliers*),
and keep the model with the largest consensus set. A final least-squares refit on those
inliers sharpens the estimate. Ideal when a clear majority-or-minority linear structure
is buried in gross contamination. Reproducible via a seeded LCG. Pure standard library.
"""


def _lcg(seed):
    state = seed & 0x7FFFFFFF
    def rand():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return state / 0x80000000
    return rand


def _ols_line(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((xs[i] - mx) ** 2 for i in range(n))
    if sxx == 0:
        return None
    sxy = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    slope = sxy / sxx
    return slope, my - slope * mx


def ransac_line(x, y, threshold, n_iterations=200, seed=1234567):
    """RANSAC line fit; returns a dict with the consensus model.

    ``threshold`` is the maximum residual for a point to count as an inlier. Runs
    ``n_iterations`` minimal (2-point) fits, keeps the model with the most inliers, and
    refits least squares on that inlier set. Returns ``slope``, ``intercept``,
    ``inliers`` (index list) and ``n_inliers``. Deterministic for a fixed ``seed``.
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have equal length")
    if n < 2:
        raise ValueError("need at least 2 points")
    if threshold <= 0:
        raise ValueError("threshold must be positive")
    rand = _lcg(seed)

    best_inliers = []
    best_model = None
    for _ in range(n_iterations):
        i = int(rand() * n) % n
        j = int(rand() * n) % n
        if i == j or x[i] == x[j]:
            continue
        slope = (y[j] - y[i]) / (x[j] - x[i])
        intercept = y[i] - slope * x[i]
        inliers = [k for k in range(n)
                   if abs(y[k] - (slope * x[k] + intercept)) <= threshold]
        if len(inliers) > len(best_inliers):
            best_inliers = inliers
            best_model = (slope, intercept)

    if best_model is None:
        raise ValueError("no valid model found (all sampled pairs share an x?)")

    # Refit least squares on the consensus set.
    if len(best_inliers) >= 2:
        xs = [x[k] for k in best_inliers]
        ys = [y[k] for k in best_inliers]
        refit = _ols_line(xs, ys)
        if refit is not None:
            best_model = refit
    slope, intercept = best_model
    return {"slope": slope, "intercept": intercept,
            "inliers": best_inliers, "n_inliers": len(best_inliers)}
