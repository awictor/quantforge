"""AAA algorithm: adaptive rational approximation (Nakatsukasa-Sete-Trefethen 2018).

AAA ("adaptive Antoulas-Anderson") is the modern method for approximating a function by a
*rational* function from sampled data. Unlike a fixed-degree fit it greedily grows a set of
support points, at each step adding the sample where the current approximant is worst, then
solving a small least-squares problem (via the SVD of a Loewner matrix) for the barycentric
weights that best fit the rest. It handles functions with poles, steep gradients, and branch
behaviour that polynomials cannot, and it is the engine behind modern rational filter design and
model-order reduction.

The result is a callable barycentric rational interpolant. Pure standard library (SVD from
:mod:`quantforge.svd`).
"""

from .svd import svd


def aaa(xs, ys, tol=1e-13, max_terms=100):
    """Rational approximation of data ``(xs, ys)`` by the AAA algorithm.

    ``xs``/``ys`` are the sample abscissae and values. Greedily selects support points until the
    max residual falls below ``tol`` (relative to ``max|ys|``) or ``max_terms`` is reached.
    Returns a callable ``r(x)`` -- the barycentric rational approximant -- together with the
    chosen support points, values and weights as attributes ``r.support_x``, ``r.support_y``,
    ``r.weights``.
    """
    m = len(xs)
    if m != len(ys):
        raise ValueError("xs and ys must have equal length")
    xs = [float(v) for v in xs]
    ys = [float(v) for v in ys]
    ymax = max(abs(y) for y in ys) or 1.0

    support = []          # indices chosen as support points
    remaining = list(range(m))
    # start with the mean as the trivial approximant
    mean_y = sum(ys) / m
    R = [mean_y] * m

    weights = []
    sx, sy = [], []
    for _ in range(min(max_terms, m)):
        # pick the remaining point with the largest residual
        j = max(remaining, key=lambda i: abs(ys[i] - R[i]))
        support.append(j)
        remaining.remove(j)
        sx = [xs[i] for i in support]
        sy = [ys[i] for i in support]

        # Loewner matrix over the remaining points
        A = []
        for i in remaining:
            row = [(ys[i] - sy[k]) / (xs[i] - sx[k]) for k in range(len(support))]
            A.append(row)
        if not A:
            weights = [1.0] * len(support)
            break
        # weights = right singular vector of smallest singular value
        _, S, V = svd(A)
        # V columns are right-singular vectors; smallest singular value is last
        n = len(support)
        smin_idx = min(range(len(S)), key=lambda t: S[t])
        weights = [V[k][smin_idx] for k in range(n)]

        # evaluate the barycentric rational at all points to get residuals
        R = _bary_eval_all(xs, sx, sy, weights)
        err = max(abs(ys[i] - R[i]) for i in range(m))
        if err <= tol * ymax:
            break

    def r(x):
        return _bary_eval_one(x, sx, sy, weights)

    r.support_x = sx
    r.support_y = sy
    r.weights = weights
    return r


def _bary_eval_one(x, sx, sy, w):
    num = 0.0
    den = 0.0
    for k in range(len(sx)):
        d = x - sx[k]
        if d == 0.0:
            return sy[k]
        t = w[k] / d
        num += t * sy[k]
        den += t
    if den == 0.0:
        return num  # degenerate; return numerator
    return num / den


def _bary_eval_all(xs, sx, sy, w):
    return [_bary_eval_one(x, sx, sy, w) for x in xs]
