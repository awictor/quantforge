"""Categorical feature encoding."""

import pytest

from quantforge import (
    fit_label_encoder, label_encode, label_decode, one_hot_encode,
)


def test_categories_sorted_and_encoded():
    enc = fit_label_encoder(["b", "a", "c", "a"])
    assert enc["categories"] == ["a", "b", "c"]
    assert label_encode(enc, ["a", "b", "c"]) == [0, 1, 2]


def test_round_trip():
    enc = fit_label_encoder(["b", "a", "c"])
    codes = label_encode(enc, ["c", "a", "b"])
    assert label_decode(enc, codes) == ["c", "a", "b"]


def test_unseen_category_maps_to_minus_one():
    enc = fit_label_encoder(["a", "b"])
    assert label_encode(enc, ["z"]) == [-1]
    assert label_decode(enc, [-1]) == [None]


def test_one_hot_rows_sum_to_one():
    enc = fit_label_encoder(["a", "b", "c"])
    oh = one_hot_encode(enc, ["a", "b", "c"])
    assert oh == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    assert all(sum(r) == 1 for r in oh)


def test_one_hot_unseen_is_zero_row():
    enc = fit_label_encoder(["a", "b"])
    assert one_hot_encode(enc, ["z"]) == [[0, 0]]


def test_validation():
    with pytest.raises(ValueError):
        fit_label_encoder([])
