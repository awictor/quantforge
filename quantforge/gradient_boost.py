"""Gradient-boosted regression trees (squared-error loss).

Boosting builds an additive model by fitting each new tree to the *residuals* of the
current ensemble, so successive shallow trees correct one another's errors -- a strong
learner assembled from weak ones. With squared-error loss the residual *is* the negative
gradient, so this is gradient descent in function space: start from the mean, then add
``learning_rate * tree(residual)`` each round. Shallow trees (``max_depth`` 2-3) plus a
small learning rate is the classic bias-variance sweet spot. Builds on
:mod:`quantforge.regression_tree`. Pure standard library.
"""

from .regression_tree import fit_regression_tree, predict_regression_tree


def fit_gradient_boost(X, y, n_estimators=100, learning_rate=0.1,
                       max_depth=3, min_samples=2):
    """Fit a gradient-boosted regression-tree ensemble (squared-error loss).

    Returns a dict with the initial ``base`` prediction (the mean of ``y``), the list
    of ``trees``, and ``learning_rate``. Each tree is fit to the current residuals and
    contributes ``learning_rate * its prediction``. Feed the result to
    :func:`predict_gradient_boost`.
    """
    n = len(y)
    if n != len(X):
        raise ValueError("X and y must have equal length")
    if n == 0:
        raise ValueError("need at least one sample")
    if n_estimators < 1:
        raise ValueError("n_estimators must be >= 1")

    base = sum(y) / n
    pred = [base] * n
    trees = []
    for _ in range(n_estimators):
        residual = [y[i] - pred[i] for i in range(n)]
        tree = fit_regression_tree(X, residual, max_depth=max_depth,
                                   min_samples=min_samples)
        update = predict_regression_tree(tree, X)
        pred = [pred[i] + learning_rate * update[i] for i in range(n)]
        trees.append(tree)
    return {"base": base, "trees": trees, "learning_rate": learning_rate}


def predict_gradient_boost(model, X_query):
    """Predict targets for ``X_query`` with a fitted gradient-boosted ensemble."""
    base = model["base"]
    lr = model["learning_rate"]
    out = [base] * len(X_query)
    for tree in model["trees"]:
        upd = predict_regression_tree(tree, X_query)
        out = [out[i] + lr * upd[i] for i in range(len(X_query))]
    return out
