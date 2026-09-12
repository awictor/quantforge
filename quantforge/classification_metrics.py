"""Binary-classification metrics.

Evaluate a probabilistic classifier's scores against binary labels:

  * ``roc_auc`` -- area under the ROC curve, computed as the Mann-Whitney rank
    statistic ``P(score of a positive > score of a negative)``. 1 = perfect,
    0.5 = random, ties count as half.
  * ``confusion_matrix`` -- (tp, fp, fn, tn) at a probability threshold.
  * ``precision_recall_f1`` -- the three summary rates at a threshold.
  * ``log_loss`` -- mean negative log-likelihood (cross-entropy); 0 = perfect.
  * ``brier_score`` -- mean squared error of the probabilities; 0 = perfect.

Pure standard library.
"""

import math


def _check(y_true, y_score):
    if len(y_true) != len(y_score):
        raise ValueError("y_true and y_score must have the same length")
    if len(y_true) == 0:
        raise ValueError("need at least one observation")
    if any(v not in (0, 1, 0.0, 1.0) for v in y_true):
        raise ValueError("y_true must be binary (0 or 1)")


def roc_auc(y_true, y_score):
    """Area under the ROC curve via the Mann-Whitney rank statistic.

    Equals the probability that a randomly chosen positive scores higher than a
    randomly chosen negative (ties = 0.5). Raises if the labels are all one class
    (AUC undefined). 1 = perfect ranking, 0.5 = random.
    """
    _check(y_true, y_score)
    pos = [y_score[i] for i in range(len(y_true)) if y_true[i] == 1]
    neg = [y_score[i] for i in range(len(y_true)) if y_true[i] == 0]
    if not pos or not neg:
        raise ValueError("AUC undefined: need both classes present")
    # Rank-sum: average rank of positives among all scores.
    order = sorted(range(len(y_score)), key=lambda i: y_score[i])
    ranks = [0.0] * len(y_score)
    i = 0
    n = len(y_score)
    while i < n:
        j = i
        while j + 1 < n and y_score[order[j + 1]] == y_score[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    rank_sum_pos = sum(ranks[i] for i in range(len(y_true)) if y_true[i] == 1)
    n_pos, n_neg = len(pos), len(neg)
    return (rank_sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def confusion_matrix(y_true, y_score, threshold=0.5):
    """Confusion counts ``(tp, fp, fn, tn)`` at a probability ``threshold``."""
    _check(y_true, y_score)
    tp = fp = fn = tn = 0
    for i in range(len(y_true)):
        pred = 1 if y_score[i] >= threshold else 0
        if y_true[i] == 1 and pred == 1:
            tp += 1
        elif y_true[i] == 0 and pred == 1:
            fp += 1
        elif y_true[i] == 1 and pred == 0:
            fn += 1
        else:
            tn += 1
    return tp, fp, fn, tn


def precision_recall_f1(y_true, y_score, threshold=0.5):
    """Precision, recall, and F1 at a threshold, as a ``(p, r, f1)`` tuple.

    Precision and recall are 0 when their denominators vanish (no predicted or no
    actual positives).
    """
    tp, fp, fn, _ = confusion_matrix(y_true, y_score, threshold)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = (2.0 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)
    return precision, recall, f1


def log_loss(y_true, y_score, eps=1e-15):
    """Mean binary cross-entropy. 0 for a perfect confident classifier."""
    _check(y_true, y_score)
    total = 0.0
    for i in range(len(y_true)):
        p = min(max(y_score[i], eps), 1.0 - eps)
        total += y_true[i] * math.log(p) + (1.0 - y_true[i]) * math.log(1.0 - p)
    return -total / len(y_true)


def brier_score(y_true, y_score):
    """Mean squared error of the predicted probabilities. 0 = perfect."""
    _check(y_true, y_score)
    return sum((y_score[i] - y_true[i]) ** 2 for i in range(len(y_true))) / len(y_true)
