"""Feature scaling: standardize, min-max, and robust scalers.

Fit a scaler's parameters on training data, then apply the same transform to test
data -- the discipline that prevents leakage. Each ``fit_*`` returns a params dict;
``transform`` and ``inverse_transform`` apply and undo it column-wise on a matrix
of rows.

  * standardize -- subtract the mean, divide by the standard deviation (z-score).
  * min-max     -- rescale each column to [0, 1].
  * robust      -- subtract the median, divide by the IQR (outlier-resistant).

Pure standard library.
"""

import math


def _columns(X):
    n = len(X)
    if n == 0:
        raise ValueError("need at least one row")
    p = len(X[0])
    if any(len(r) != p for r in X):
        raise ValueError("all rows must have the same length")
    return n, p


def _std(col):
    n = len(col)
    m = sum(col) / n
    return math.sqrt(sum((v - m) ** 2 for v in col) / n)


def _quantile(sorted_col, q):
    n = len(sorted_col)
    if n == 1:
        return sorted_col[0]
    pos = q * (n - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 < n:
        return sorted_col[lo] * (1.0 - frac) + sorted_col[lo + 1] * frac
    return sorted_col[lo]


def fit_standardize(X):
    """Fit a z-score scaler: per-column ``(mean, std)`` (zero std -> 1)."""
    n, p = _columns(X)
    centers = [sum(X[t][j] for t in range(n)) / n for j in range(p)]
    scales = [_std([X[t][j] for t in range(n)]) or 1.0 for j in range(p)]
    return {"kind": "standardize", "center": centers, "scale": scales}


def fit_min_max(X):
    """Fit a min-max scaler: per-column ``(min, range)`` (zero range -> 1)."""
    n, p = _columns(X)
    mins = [min(X[t][j] for t in range(n)) for j in range(p)]
    maxs = [max(X[t][j] for t in range(n)) for j in range(p)]
    scales = [(maxs[j] - mins[j]) or 1.0 for j in range(p)]
    return {"kind": "min_max", "center": mins, "scale": scales}


def fit_robust(X):
    """Fit a robust scaler: per-column ``(median, IQR)`` (zero IQR -> 1)."""
    n, p = _columns(X)
    centers = []
    scales = []
    for j in range(p):
        col = sorted(X[t][j] for t in range(n))
        centers.append(_quantile(col, 0.5))
        iqr = _quantile(col, 0.75) - _quantile(col, 0.25)
        scales.append(iqr or 1.0)
    return {"kind": "robust", "center": centers, "scale": scales}


def transform(params, X):
    """Apply a fitted scaler: ``(x - center) / scale`` column-wise."""
    n, p = _columns(X)
    c, s = params["center"], params["scale"]
    if len(c) != p:
        raise ValueError("params do not match the number of columns")
    return [[(X[t][j] - c[j]) / s[j] for j in range(p)] for t in range(n)]


def inverse_transform(params, X):
    """Undo a fitted scaler: ``x * scale + center`` column-wise."""
    n, p = _columns(X)
    c, s = params["center"], params["scale"]
    if len(c) != p:
        raise ValueError("params do not match the number of columns")
    return [[X[t][j] * s[j] + c[j] for j in range(p)] for t in range(n)]
