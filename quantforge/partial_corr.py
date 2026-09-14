"""Partial and semi-partial correlation.

The correlation between two variables after removing the linear influence of one or more
control variables -- the tool for asking "do X and Y move together *once we account for*
Z?" A spurious correlation driven entirely by a common cause drops to zero when that
cause is controlled. Computed by regressing out the controls (ordinary least squares) and
correlating the residuals. Pure standard library.
"""

import math


def _mean(v):
    return sum(v) / len(v)


def _pearson(a, b):
    n = len(a)
    ma, mb = _mean(a), _mean(b)
    sab = sum((a[i] - ma) * (b[i] - mb) for i in range(n))
    saa = sum((a[i] - ma) ** 2 for i in range(n))
    sbb = sum((b[i] - mb) ** 2 for i in range(n))
    if saa <= 0.0 or sbb <= 0.0:
        raise ValueError("a variable has zero variance")
    return sab / math.sqrt(saa * sbb)


def _ols_residuals(y, controls):
    """Residuals of ``y`` after OLS on the control columns (with an intercept).

    ``controls`` is a list of predictor columns (each a list). Solves the normal
    equations for the intercept + control coefficients and returns ``y - y_hat``.
    """
    n = len(y)
    # Design matrix rows: [1, c1_i, c2_i, ...].
    p = len(controls)
    X = [[1.0] + [controls[j][i] for j in range(p)] for i in range(n)]
    m = p + 1
    # Normal equations (X'X) beta = X'y.
    xtx = [[sum(X[k][a] * X[k][b] for k in range(n)) for b in range(m)] for a in range(m)]
    xty = [sum(X[k][a] * y[k] for k in range(n)) for a in range(m)]
    beta = _solve(xtx, xty)
    resid = [y[i] - sum(beta[a] * X[i][a] for a in range(m)) for i in range(n)]
    return resid


def _solve(A, b):
    """Gaussian elimination with partial pivoting for a small dense system."""
    n = len(A)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-15:
            raise ValueError("singular control matrix (collinear controls?)")
        M[col], M[piv] = M[piv], M[col]
        pivot = M[col][col]
        for r in range(n):
            if r != col:
                f = M[r][col] / pivot
                for c in range(col, n + 1):
                    M[r][c] -= f * M[col][c]
    return [M[i][n] / M[i][i] for i in range(n)]


def partial_correlation(x, y, controls):
    """Partial correlation of ``x`` and ``y`` controlling for ``controls``.

    ``controls`` is a single control column (list) or a list of control columns. Removes
    the linear effect of the controls from both ``x`` and ``y`` (OLS) and correlates the
    residuals. Returns a coefficient in ``[-1, 1]``; a correlation that is purely due to
    the controls drops toward 0.
    """
    if controls and not isinstance(controls[0], (list, tuple)):
        controls = [controls]
    cols = [list(c) for c in controls]
    rx = _ols_residuals(list(x), cols)
    ry = _ols_residuals(list(y), cols)
    return _pearson(rx, ry)


def semipartial_correlation(x, y, controls):
    """Semi-partial (part) correlation: control the ``controls`` out of ``y`` only.

    Correlates raw ``x`` with the residual of ``y`` after regressing out the controls --
    the unique contribution of ``x`` to ``y`` beyond the controls. Returns a coefficient
    in ``[-1, 1]``.
    """
    if controls and not isinstance(controls[0], (list, tuple)):
        controls = [controls]
    cols = [list(c) for c in controls]
    ry = _ols_residuals(list(y), cols)
    return _pearson(list(x), ry)
