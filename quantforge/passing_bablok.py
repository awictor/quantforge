"""Passing-Bablok regression: robust nonparametric method comparison.

The standard method-comparison fit in clinical chemistry: like Deming it treats both
axes as error-prone, but distribution-free and robust to outliers. It takes the median
of all pairwise slopes ``(y_j - y_i)/(x_j - x_i)`` -- but with two corrections that make
it consistent for method comparison rather than a plain Theil-Sen fit:

  * slopes equal to ``-1`` (and undefined ones from tied ``x``) are discarded, and
  * the median is *shifted* by ``K``, the number of slopes below ``-1``, which removes
    the bias the discarded/negative slopes would otherwise introduce.

The intercept is then ``median(y_i - b x_i)``. Pure standard library; ``O(n^2)``.
"""


def _median_sorted(values):
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2:
        return s[mid]
    return 0.5 * (s[mid - 1] + s[mid])


def passing_bablok_regression(x, y):
    """Passing-Bablok regression of ``y`` on ``x``.

    Returns ``(slope, intercept)`` for ``y = slope * x + intercept``. Robust to
    outliers in either variable and symmetric in ``x`` and ``y`` (up to reciprocal
    slope), the standard nonparametric alternative to Deming regression.
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have equal length")
    if n < 2:
        raise ValueError("need at least 2 points")

    slopes = []
    k_below = 0             # count of slopes strictly less than -1
    for i in range(n):
        for j in range(i + 1, n):
            dx = x[j] - x[i]
            dy = y[j] - y[i]
            if dx == 0 and dy == 0:
                continue                    # identical points carry no slope
            if dx == 0:
                continue                    # vertical pair: undefined slope, skip
            s = dy / dx
            if s == -1.0:
                continue                    # excluded by the Passing-Bablok rule
            slopes.append(s)
            if s < -1.0:
                k_below += 1

    if not slopes:
        raise ValueError("no valid pairwise slopes")

    slopes.sort()
    n_s = len(slopes)

    def at(idx):
        # Clamp into range: the K shift assumes a positive (method-comparison)
        # slope, so idx stays in range there; clamping keeps steep-negative data
        # from indexing past the end rather than crashing.
        if idx < 0:
            idx = 0
        elif idx >= n_s:
            idx = n_s - 1
        return slopes[idx]

    # Shifted median: offset the rank by the number of slopes below -1.
    if n_s % 2 == 1:
        slope = at((n_s - 1) // 2 + k_below)
    else:
        slope = 0.5 * (at(n_s // 2 - 1 + k_below) + at(n_s // 2 + k_below))

    intercept = _median_sorted([y[i] - slope * x[i] for i in range(n)])
    return slope, intercept
