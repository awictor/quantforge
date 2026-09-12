"""Cross-validation splitters and scoring."""

import random

import pytest

from quantforge import k_fold_indices, train_test_split, cross_val_score


def test_folds_partition_exactly():
    folds = list(k_fold_indices(23, 5))
    all_test = [i for _, te in folds for i in te]
    assert sorted(all_test) == list(range(23))
    assert len(all_test) == len(set(all_test))


def test_train_test_complementary_and_disjoint():
    for tr, te in k_fold_indices(23, 5):
        assert sorted(tr + te) == list(range(23))
        assert set(tr).isdisjoint(te)


def test_fold_sizes_balanced():
    sizes = [len(te) for _, te in k_fold_indices(23, 5)]
    assert max(sizes) - min(sizes) <= 1
    assert len(sizes) == 5


def test_shuffle_reproducible_and_seed_sensitive():
    a = [te for _, te in k_fold_indices(20, 4, True, 42)]
    b = [te for _, te in k_fold_indices(20, 4, True, 42)]
    c = [te for _, te in k_fold_indices(20, 4, True, 99)]
    assert a == b
    assert a != c


def test_train_test_split_sizes():
    tr, te = train_test_split(100, 0.2)
    assert len(tr) == 80 and len(te) == 20
    assert sorted(tr + te) == list(range(100))


def test_cross_val_score_runs_over_folds():
    rng = random.Random(1)
    X = [[i] for i in range(50)]
    y = [rng.gauss(0, 1) for _ in range(50)]

    def fit(Xtr, ytr):
        return sum(ytr) / len(ytr)          # constant mean model

    def score(m, Xte, yte):
        return sum((m - v) ** 2 for v in yte) / len(yte)

    scores = cross_val_score(X, y, fit, score, 5)
    assert len(scores) == 5
    assert all(s >= 0 for s in scores)


def test_validation():
    with pytest.raises(ValueError):
        list(k_fold_indices(3, 5))          # n < k
    with pytest.raises(ValueError):
        list(k_fold_indices(10, 1))         # k < 2
    with pytest.raises(ValueError):
        train_test_split(10, 1.5)           # fraction out of range
