"""CART classification tree (recursive Gini splitting).

Extends the decision stump to a full tree: recursively split each node on the
feature/threshold that most reduces Gini impurity, stopping at a maximum depth, a
minimum node size, or a pure node. Each leaf predicts its majority class. Greedy
CART -- the standard interpretable classifier and the building block of random
forests. Pure standard library.
"""


def _gini(labels):
    n = len(labels)
    counts = {}
    for y in labels:
        counts[y] = counts.get(y, 0) + 1
    return 1.0 - sum((c / n) ** 2 for c in counts.values())


def _majority(labels):
    counts = {}
    for y in labels:
        counts[y] = counts.get(y, 0) + 1
    return max(counts, key=lambda k: counts[k])


def _best_split(X, y, idx):
    n = len(idx)
    d = len(X[0])
    best = None
    best_g = float("inf")
    for f in range(d):
        vals = sorted(set(X[i][f] for i in idx))
        for a, b in zip(vals, vals[1:]):
            thr = 0.5 * (a + b)
            left = [i for i in idx if X[i][f] <= thr]
            right = [i for i in idx if X[i][f] > thr]
            if not left or not right:
                continue
            g = (len(left) * _gini([y[i] for i in left])
                 + len(right) * _gini([y[i] for i in right])) / n
            if g < best_g:
                best_g = g
                best = (f, thr, left, right)
    return best


def _build(X, y, idx, depth, max_depth, min_samples):
    labels = [y[i] for i in idx]
    node = {"prediction": _majority(labels)}
    if depth >= max_depth or len(idx) < min_samples or _gini(labels) == 0.0:
        return node
    split = _best_split(X, y, idx)
    if split is None:
        return node
    f, thr, left, right = split
    node["feature"] = f
    node["threshold"] = thr
    node["left"] = _build(X, y, left, depth + 1, max_depth, min_samples)
    node["right"] = _build(X, y, right, depth + 1, max_depth, min_samples)
    return node


def fit_decision_tree(X, y, max_depth=5, min_samples=2):
    """Fit a CART classification tree.

    Parameters
    ----------
    X : list[list[float]]
        ``n`` rows of ``d`` features.
    y : list
        Class labels.
    max_depth : int
        Maximum tree depth (>= 1). Depth 1 is a decision stump.
    min_samples : int
        Minimum samples required to attempt a split.

    Returns
    -------
    dict
        Nested tree; internal nodes have ``feature``/``threshold``/``left``/
        ``right``, leaves have only ``prediction``.
    """
    n = len(X)
    if n == 0:
        raise ValueError("need at least one observation")
    if n != len(y):
        raise ValueError("X and y must have the same length")
    if max_depth < 1:
        raise ValueError("max_depth must be >= 1")
    return _build(X, y, list(range(n)), 0, max_depth, min_samples)


def predict_decision_tree(tree, X_query):
    """Predict labels for query rows by walking the tree to a leaf."""
    def one(node, x):
        while "feature" in node:
            node = node["left"] if x[node["feature"]] <= node["threshold"] else node["right"]
        return node["prediction"]
    return [one(tree, x) for x in X_query]


def tree_depth(tree):
    """Depth of a fitted tree (a single leaf has depth 0)."""
    if "feature" not in tree:
        return 0
    return 1 + max(tree_depth(tree["left"]), tree_depth(tree["right"]))
