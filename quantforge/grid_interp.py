"""Interpolation on a 2-D regular grid: bilinear and nearest-neighbour.

Given a table of values ``z[i][j]`` sampled on a rectilinear grid with coordinate axes
``xs`` (columns) and ``ys`` (rows), estimate the value at an arbitrary ``(x, y)``.
Bilinear interpolation blends the four surrounding grid values; nearest-neighbour snaps to
the closest node. Both clamp to the grid edges for out-of-range queries. This is the tool
behind lookup tables -- vol surfaces, response grids, heightmaps. Pure standard library.
"""

import bisect


def _locate(axis, v):
    """Index ``i`` and fraction ``f`` so ``v`` lies between ``axis[i]`` and ``axis[i+1]``.

    Clamps to the ends when ``v`` is outside the axis range.
    """
    n = len(axis)
    if v <= axis[0]:
        return 0, 0.0
    if v >= axis[-1]:
        return n - 2, 1.0
    j = bisect.bisect_right(axis, v) - 1
    if j >= n - 1:
        j = n - 2
    span = axis[j + 1] - axis[j]
    f = 0.0 if span == 0 else (v - axis[j]) / span
    return j, f


def bilinear_interp(xs, ys, z, x, y):
    """Bilinear interpolation of grid ``z`` at ``(x, y)``.

    ``xs`` are the column coordinates (length ``ncols``), ``ys`` the row coordinates
    (length ``nrows``), and ``z`` an ``nrows x ncols`` value table (``z[row][col]``). Both
    axes must be strictly increasing. Out-of-range queries clamp to the nearest edge.
    Exact at grid nodes; reduces to linear interpolation along a grid line.
    """
    if len(xs) < 2 or len(ys) < 2:
        raise ValueError("each axis needs at least two points")
    if len(z) != len(ys) or any(len(row) != len(xs) for row in z):
        raise ValueError("z must be shaped (len(ys), len(xs))")
    i, fx = _locate(xs, x)
    j, fy = _locate(ys, y)
    z00 = z[j][i]
    z01 = z[j][i + 1]
    z10 = z[j + 1][i]
    z11 = z[j + 1][i + 1]
    top = z00 * (1 - fx) + z01 * fx
    bot = z10 * (1 - fx) + z11 * fx
    return top * (1 - fy) + bot * fy


def nearest_interp(xs, ys, z, x, y):
    """Nearest-neighbour lookup on the grid: the value at the closest node to ``(x, y)``.

    Snaps to the nearest grid coordinate on each axis (clamping out-of-range). Piecewise
    constant, so it preserves the exact sampled values -- useful for categorical or
    quantized grids.
    """
    if len(xs) < 1 or len(ys) < 1:
        raise ValueError("each axis needs at least one point")
    if len(z) != len(ys) or any(len(row) != len(xs) for row in z):
        raise ValueError("z must be shaped (len(ys), len(xs))")
    i, fx = _locate(xs, x) if len(xs) >= 2 else (0, 0.0)
    j, fy = _locate(ys, y) if len(ys) >= 2 else (0, 0.0)
    ci = i + (1 if fx >= 0.5 else 0)
    cj = j + (1 if fy >= 0.5 else 0)
    ci = min(ci, len(xs) - 1)
    cj = min(cj, len(ys) - 1)
    return z[cj][ci]
