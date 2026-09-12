"""Polynomial and interaction feature expansion."""

from math import comb

import pytest

from quantforge import polynomial_features


def test_degree_one_with_bias_is_identity():
    F, terms = polynomial_features([[2, 3], [4, 5]], degree=1)
    assert F[0] == [1.0, 2.0, 3.0]
    assert terms == [(), (0,), (1,)]


def test_degree_two_includes_squares_and_products():
    F, _ = polynomial_features([[2, 3]], degree=2)
    # [1, x1, x2, x1^2, x1*x2, x2^2]
    assert F[0] == [1.0, 2.0, 3.0, 4.0, 6.0, 9.0]


def test_feature_count_matches_formula():
    p, d = 2, 3
    _, terms = polynomial_features([[1, 1]], d)
    assert len(terms) == comb(p + d, d)


def test_interaction_only_drops_pure_powers():
    F, _ = polynomial_features([[2, 3]], degree=2, interaction_only=True)
    assert F[0] == [1.0, 2.0, 3.0, 6.0]      # no x1^2 or x2^2


def test_no_bias():
    F, terms = polynomial_features([[2, 3]], degree=1, include_bias=False)
    assert F[0] == [2.0, 3.0]
    assert () not in terms


def test_validation():
    with pytest.raises(ValueError):
        polynomial_features([], 2)
    with pytest.raises(ValueError):
        polynomial_features([[1, 2]], 0)
    with pytest.raises(ValueError):
        polynomial_features([[1, 2], [3]], 2)   # ragged rows
