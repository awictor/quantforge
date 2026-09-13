"""Platt scaling: sigmoid probability calibration of classifier scores.

Platt (1999) calibrates a classifier's raw scores ``s`` into probabilities by fitting
a one-parameter-pair sigmoid

    P(y = 1 | s) = 1 / (1 + exp(-(A s + B))),

choosing ``A, B`` to minimize the regularized logistic loss on the training scores.
(Platt's paper writes the sigmoid with ``+(A s + B)`` in the exponent and so reports
``A < 0`` for a good classifier; the logistic-standard sign used here gives ``A > 0``
for the same fit -- the calibrated probabilities are identical either way.)
Platt's target smoothing replaces the hard 0/1 labels with ``1/(N_neg + 2)`` and
``(N_pos + 1)/(N_pos + 2)`` to avoid overfitting the tails. Unlike isotonic
calibration this assumes a sigmoidal miscalibration shape, so it needs far less data
and never overfits to a step function. Fit by Newton's method; pure standard library.
"""

import math


def _sigmoid(z):
    # Numerically stable logistic.
    if z >= 0.0:
        e = math.exp(-z)
        return 1.0 / (1.0 + e)
    e = math.exp(z)
    return e / (1.0 + e)


def platt_fit(scores, labels, max_iter=100, tol=1e-10):
    """Fit Platt sigmoid parameters ``(A, B)`` to scores and binary labels.

    ``P(y=1|s) = 1 / (1 + exp(-(A s + B)))``. Uses Platt's smoothed targets and
    Newton's method on the regularized logistic loss. Returns ``(A, B)``; a
    well-separated classifier gives ``A > 0`` (probability rises with the score) in
    this logistic-standard sign convention. Pure standard library.
    """
    n = len(scores)
    if n != len(labels):
        raise ValueError("scores and labels must have equal length")
    if n == 0:
        raise ValueError("need at least one point")
    if any(l not in (0, 1, 0.0, 1.0) for l in labels):
        raise ValueError("labels must be 0 or 1")

    prior1 = sum(1 for l in labels if l == 1)
    prior0 = n - prior1
    # Platt's smoothed targets.
    hi = (prior1 + 1.0) / (prior1 + 2.0)
    lo = 1.0 / (prior0 + 2.0)
    t = [hi if l == 1 else lo for l in labels]

    A = 0.0
    B = math.log((prior0 + 1.0) / (prior1 + 1.0))
    lam = 1e-12
    for _ in range(max_iter):
        # Gradient and Hessian of the log loss w.r.t. (A, B).
        g1 = g2 = h11 = h12 = h22 = 0.0
        for i in range(n):
            z = A * scores[i] + B
            p = _sigmoid(z)          # p = P(y=1) = 1/(1+exp(-z))
            # NLL gradient w.r.t. z is (p - t); Newton solves H dx = -grad, and g1/g2
            # below accumulate -grad = (t - p) times [s, 1].
            d1 = t[i] - p            # (target - model prob)
            wp = p * (1.0 - p)
            g1 += scores[i] * d1
            g2 += d1
            h11 += scores[i] * scores[i] * wp
            h12 += scores[i] * wp
            h22 += wp
        # Newton step solving H dx = g (with tiny ridge for stability).
        h11 += lam
        h22 += lam
        det = h11 * h22 - h12 * h12
        if abs(det) < 1e-300:
            break
        dA = (h22 * g1 - h12 * g2) / det
        dB = (h11 * g2 - h12 * g1) / det
        A += dA
        B += dB
        if abs(dA) < tol and abs(dB) < tol:
            break
    return A, B


def platt_predict(scores, A, B):
    """Apply a fitted Platt sigmoid to scores, returning calibrated probabilities."""
    return [_sigmoid(A * s + B) for s in scores]


def platt_calibrate(train_scores, train_labels, new_scores=None, max_iter=100):
    """Fit Platt scaling and return ``(params, predict)``.

    ``params`` is ``(A, B)`` and ``predict(scores)`` maps raw scores to calibrated
    probabilities. If ``new_scores`` is given, returns ``(params, calibrated_new)``
    instead, applying the fit directly.
    """
    A, B = platt_fit(train_scores, train_labels, max_iter=max_iter)
    if new_scores is not None:
        return (A, B), platt_predict(new_scores, A, B)

    def predict(scores):
        return platt_predict(scores, A, B)

    return (A, B), predict
