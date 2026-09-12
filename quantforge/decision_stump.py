"""Decision stump: a one-split classifier via Gini impurity.

A stump is a depth-1 decision tree -- it picks the single feature and threshold
whose split most reduces class impurity, then predicts the majority class on each
side. Gini impurity of a node is ``1 - sum_c p_c^2`` (0 = pure). The best split
minimizes the size-weighted Gini of the two children over every feature and every
midpoint between sorted distinct values. Simple, interpretable, and the base
learner of boosting. Pure standard library.
"""


def gini_impurity(labels):
    """Gini impurity ``1 - sum_c p_c^2`` of a label list (0 = pure)."""
    n = len(labels)
    if n == 0:
        raise ValueError("need at least one label")
    counts = {}
    for y in labels:
        counts[y] = counts.get(y, 0) + 1
    return 1.0 - sum((c / n) ** 2 for c in counts.values())


def _majority(labels):
    counts = {}
    for y in labels:
        counts[y] = counts.get(y, 0) + 1
    return max(counts, key=lambda k: counts[k])


def fit_decision_stump(X, y):
    """Fit a decision stump: the best single-feature threshold split.

    Scans every feature and every midpoint between adjacent sorted values, scoring
    each candidate by the size-weighted Gini of the two sides. Returns
    ``{"feature", "threshold", "left_label", "right_label", "gini"}`` where points
    with ``x[feature] <= threshold`` take ``left_label``.
    """
    n = len(X)
    if n == 0:
        raise ValueError("need at least one observation")
    if n != len(y):
        raise ValueError("X and y must have the same length")
    d = len(X[0])

    best = None
    best_gini = float("inf")
    for f in range(d):
        vals = sorted(set(X[i][f] for i in range(n)))
        for a, b in zip(vals, vals[1:]):
            thr = 0.5 * (a + b)
            left = [y[i] for i in range(n) if X[i][f] <= thr]
            right = [y[i] for i in range(n) if X[i][f] > thr]
            if not left or not right:
                continue
            w = (len(left) * gini_impurity(left)
                 + len(right) * gini_impurity(right)) / n
            if w < best_gini:
                best_gini = w
                best = {
                    "feature": f,
                    "threshold": thr,
                    "left_label": _majority(left),
                    "right_label": _majority(right),
                    "gini": w,
                }
    if best is None:
        # No split possible (all rows identical); predict the global majority.
        maj = _majority(y)
        best = {"feature": 0, "threshold": float("inf"),
                "left_label": maj, "right_label": maj,
                "gini": gini_impurity(y)}
    return best


def predict_decision_stump(stump, X_query):
    """Predict labels for query rows under a fitted stump."""
    f, thr = stump["feature"], stump["threshold"]
    return [stump["left_label"] if x[f] <= thr else stump["right_label"]
            for x in X_query]
