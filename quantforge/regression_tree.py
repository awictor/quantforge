"""CART regression tree (recursive variance-reduction splitting).

The regression analogue of a classification tree: instead of Gini impurity it splits to
minimize the sum of squared errors within each child, and a leaf predicts the *mean* of
the responses that fall in it. A single tree is a piecewise-constant regressor -- it
captures nonlinearities and interactions a linear model can't, at the cost of a step
function -- and is the building block of gradient-boosted and bagged regressors. Pure
standard library.
"""


def _sse(ys):
    n = len(ys)
    if n == 0:
        return 0.0
    m = sum(ys) / n
    return sum((y - m) ** 2 for y in ys)


def _best_split(X, y):
    n = len(y)
    p = len(X[0])
    parent = _sse(y)
    best = None
    best_gain = 0.0
    for j in range(p):
        # Sort by feature j.
        order = sorted(range(n), key=lambda i: X[i][j])
        xs = [X[i][j] for i in order]
        ys = [y[i] for i in order]
        for s in range(1, n):
            if xs[s] == xs[s - 1]:
                continue
            left, right = ys[:s], ys[s:]
            gain = parent - (_sse(left) + _sse(right))
            if gain > best_gain:
                best_gain = gain
                thr = 0.5 * (xs[s] + xs[s - 1])
                best = (j, thr)
    return best, best_gain


def fit_regression_tree(X, y, max_depth=5, min_samples=2):
    """Fit a CART regression tree.

    ``X`` is a list of feature rows, ``y`` the numeric targets. Splits greedily to
    maximize squared-error reduction until ``max_depth`` or ``min_samples`` stops it;
    leaves store the mean target. Returns a nested-dict tree for
    :func:`predict_regression_tree`.
    """
    n = len(y)
    if n != len(X):
        raise ValueError("X and y must have equal length")
    if n == 0:
        raise ValueError("need at least one sample")

    def build(idx, depth):
        ys = [y[i] for i in idx]
        leaf = {"leaf": True, "value": sum(ys) / len(ys)}
        if depth >= max_depth or len(idx) < 2 * min_samples or _sse(ys) == 0.0:
            return leaf
        sub_X = [X[i] for i in idx]
        split, gain = _best_split(sub_X, ys)
        if split is None or gain <= 0.0:
            return leaf
        j, thr = split
        left = [i for i in idx if X[i][j] <= thr]
        right = [i for i in idx if X[i][j] > thr]
        if len(left) < min_samples or len(right) < min_samples:
            return leaf
        return {"leaf": False, "feature": j, "threshold": thr,
                "left": build(left, depth + 1), "right": build(right, depth + 1)}

    return build(list(range(n)), 0)


def predict_regression_tree(tree, X_query):
    """Predict targets for rows ``X_query`` with a fitted regression tree."""
    def one(row):
        node = tree
        while not node["leaf"]:
            if row[node["feature"]] <= node["threshold"]:
                node = node["left"]
            else:
                node = node["right"]
        return node["value"]
    return [one(row) for row in X_query]
