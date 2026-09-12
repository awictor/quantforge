"""Random forest classifier (bagged CART trees).

Trains an ensemble of CART trees, each on a bootstrap resample of the data, and
classifies by majority vote across the trees. Bagging decorrelates the trees and
averages out the high variance of a single deep tree, so the forest generalizes
better than any one member. A deterministic linear-congruential stream makes the
bootstrap reproducible per seed. Pure standard library (builds on
``decision_tree``).
"""

from .decision_tree import fit_decision_tree, predict_decision_tree


def _lcg(seed):
    state = seed & 0x7FFFFFFF

    def _next():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return state
    return _next


def fit_random_forest(X, y, n_trees=10, max_depth=5, min_samples=2, seed=1234567):
    """Fit a random forest of bootstrap-resampled CART trees.

    Parameters
    ----------
    X, y : data and labels.
    n_trees : int
        Number of trees in the ensemble (>= 1).
    max_depth, min_samples : passed to each tree.
    seed : int
        Seed for the reproducible bootstrap resampling.

    Returns
    -------
    dict
        ``{"trees": [...], "classes": sorted labels}``.
    """
    n = len(X)
    if n == 0:
        raise ValueError("need at least one observation")
    if n != len(y):
        raise ValueError("X and y must have the same length")
    if n_trees < 1:
        raise ValueError("n_trees must be >= 1")
    rnd = _lcg(seed)
    trees = []
    for _ in range(n_trees):
        idx = [rnd() % n for _ in range(n)]         # bootstrap sample
        Xb = [X[i] for i in idx]
        yb = [y[i] for i in idx]
        trees.append(fit_decision_tree(Xb, yb, max_depth=max_depth,
                                       min_samples=min_samples))
    return {"trees": trees, "classes": sorted(set(y))}


def predict_random_forest(forest, X_query):
    """Majority-vote prediction across the forest's trees."""
    trees = forest["trees"]
    per_tree = [predict_decision_tree(t, X_query) for t in trees]
    preds = []
    for i in range(len(X_query)):
        votes = {}
        for t in range(len(trees)):
            lab = per_tree[t][i]
            votes[lab] = votes.get(lab, 0) + 1
        preds.append(max(votes, key=lambda k: votes[k]))
    return preds
