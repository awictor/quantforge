"""k-nearest-neighbor classification and regression.

A lazy, non-parametric learner: to predict a query, find its ``k`` closest
training points by Euclidean distance and either vote on their labels
(classification) or average their targets (regression). No training step -- the
model is the data. Simple, strong baseline for both tasks. Pure standard library.
"""

import math


def _dist2(a, b):
    return sum((a[i] - b[i]) ** 2 for i in range(len(a)))


def _neighbors(X_train, query, k):
    dists = sorted(range(len(X_train)), key=lambda i: _dist2(X_train[i], query))
    return dists[:k]


def knn_classify(X_train, y_train, X_query, k=3):
    """Classify each query row by majority vote of its ``k`` nearest neighbors.

    Ties are broken toward the label that is closest on average (the first
    encountered at minimum total distance). Returns a predicted label per query.
    """
    n = len(X_train)
    if n != len(y_train):
        raise ValueError("X_train and y_train must have the same length")
    if k < 1 or k > n:
        raise ValueError("k must satisfy 1 <= k <= len(X_train)")
    preds = []
    for q in X_query:
        idx = _neighbors(X_train, q, k)
        votes = {}
        dsum = {}
        for i in idx:
            lab = y_train[i]
            votes[lab] = votes.get(lab, 0) + 1
            dsum[lab] = dsum.get(lab, 0.0) + math.sqrt(_dist2(X_train[i], q))
        # Highest vote, ties broken by smallest total distance.
        best = min(votes, key=lambda lab: (-votes[lab], dsum[lab]))
        preds.append(best)
    return preds


def knn_regress(X_train, y_train, X_query, k=3):
    """Predict each query as the mean target of its ``k`` nearest neighbors."""
    n = len(X_train)
    if n != len(y_train):
        raise ValueError("X_train and y_train must have the same length")
    if k < 1 or k > n:
        raise ValueError("k must satisfy 1 <= k <= len(X_train)")
    preds = []
    for q in X_query:
        idx = _neighbors(X_train, q, k)
        preds.append(sum(y_train[i] for i in idx) / k)
    return preds
