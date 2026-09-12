"""Binary-classification metrics."""

import random

import pytest

from quantforge import (
    roc_auc, confusion_matrix, precision_recall_f1, log_loss, brier_score,
)


def test_perfect_and_reversed_auc():
    yt = [0, 0, 0, 1, 1, 1]
    ys = [0.1, 0.2, 0.3, 0.7, 0.8, 0.9]
    assert roc_auc(yt, ys) == 1.0
    assert roc_auc(yt, [0.9, 0.8, 0.7, 0.3, 0.2, 0.1]) == 0.0


def test_random_auc_near_half():
    rng = random.Random(1)
    y = [rng.randint(0, 1) for _ in range(5000)]
    s = [rng.random() for _ in range(5000)]
    assert abs(roc_auc(y, s) - 0.5) < 0.03


def test_ties_count_half():
    assert roc_auc([0, 1], [0.5, 0.5]) == 0.5


def test_confusion_and_prf_perfect():
    yt = [0, 0, 0, 1, 1, 1]
    ys = [0.1, 0.2, 0.3, 0.7, 0.8, 0.9]
    assert confusion_matrix(yt, ys, 0.5) == (3, 0, 0, 3)
    p, r, f = precision_recall_f1(yt, ys)
    assert p == 1.0 and r == 1.0 and f == 1.0


def test_log_loss_and_brier():
    yt = [0, 0, 0, 1, 1, 1]
    assert log_loss(yt, [0.01, 0.02, 0.03, 0.97, 0.98, 0.99]) < 0.05
    assert brier_score([1, 0], [1.0, 0.0]) == 0.0
    assert abs(brier_score([1, 0, 1, 0], [0.5] * 4) - 0.25) < 1e-12


def test_validation():
    with pytest.raises(ValueError):
        roc_auc([1, 1, 1], [0.5, 0.6, 0.7])       # one class only
    with pytest.raises(ValueError):
        log_loss([0.5], [0.5])                    # non-binary label
    with pytest.raises(ValueError):
        brier_score([1], [0.5, 0.6])              # length mismatch
