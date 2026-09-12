"""Gaussian naive Bayes classifier.

Assumes the features are conditionally independent given the class and normally
distributed within each class. Training estimates each class's prior and, per
feature, its mean and variance. Prediction picks the class maximizing the log
posterior

    log P(c) + sum_j log N(x_j; mu_{c,j}, var_{c,j}),

computed in log space for numerical stability. Despite the crude independence
assumption it is a fast, strong baseline. Pure standard library.
"""

import math


def fit_gaussian_nb(X, y, var_smoothing=1e-9):
    """Fit a Gaussian naive Bayes model.

    Parameters
    ----------
    X : list[list[float]]
        ``n`` rows of ``d`` features.
    y : list
        Class labels (any hashable).
    var_smoothing : float
        Added to every variance for numerical stability (avoids zero variance on
        constant features).

    Returns
    -------
    dict
        ``classes`` (sorted), ``priors``, ``means`` and ``variances`` (per class,
        per feature), keyed by class label.
    """
    n = len(X)
    if n == 0:
        raise ValueError("need at least one observation")
    if n != len(y):
        raise ValueError("X and y must have the same length")
    d = len(X[0])
    classes = sorted(set(y))
    by_class = {c: [X[i] for i in range(n) if y[i] == c] for c in classes}

    priors, means, variances = {}, {}, {}
    for c in classes:
        rows = by_class[c]
        m = len(rows)
        priors[c] = m / n
        means[c] = [sum(r[j] for r in rows) / m for j in range(d)]
        variances[c] = [
            sum((r[j] - means[c][j]) ** 2 for r in rows) / m + var_smoothing
            for j in range(d)
        ]
    return {"classes": classes, "priors": priors, "means": means,
            "variances": variances}


def _log_likelihood(model, x, c):
    ll = math.log(model["priors"][c])
    mu = model["means"][c]
    var = model["variances"][c]
    for j in range(len(x)):
        ll += -0.5 * (math.log(2.0 * math.pi * var[j])
                      + (x[j] - mu[j]) ** 2 / var[j])
    return ll


def predict_gaussian_nb(model, X_query):
    """Predict the most probable class for each query row."""
    preds = []
    for x in X_query:
        best = max(model["classes"], key=lambda c: _log_likelihood(model, x, c))
        preds.append(best)
    return preds


def predict_proba_gaussian_nb(model, X_query):
    """Posterior class probabilities per query (softmax of the log posteriors)."""
    out = []
    for x in X_query:
        lls = [_log_likelihood(model, x, c) for c in model["classes"]]
        mx = max(lls)
        exps = [math.exp(v - mx) for v in lls]
        z = sum(exps)
        out.append({model["classes"][i]: exps[i] / z
                    for i in range(len(model["classes"]))})
    return out
