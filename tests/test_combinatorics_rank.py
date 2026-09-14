"""Combinatorial ranking: Gray code, permutation and combination rank/unrank."""

import itertools
from math import comb

import pytest

from quantforge import (
    gray_code,
    gray_decode,
    permutation_unrank,
    permutation_rank,
    combination_unrank,
    combination_rank,
)


def test_gray_code():
    for n in range(1, 2000):
        assert gray_decode(gray_code(n)) == n
        diff = gray_code(n) ^ gray_code(n - 1)
        assert diff != 0 and diff & (diff - 1) == 0     # exactly one bit differs
    assert [gray_code(i) for i in range(8)] == [0, 1, 3, 2, 6, 7, 5, 4]


def test_permutation_rank_unrank():
    for n in range(1, 7):
        perms = list(itertools.permutations(range(n)))
        for r, p in enumerate(perms):
            assert tuple(permutation_unrank(r, n)) == p
            assert permutation_rank(list(p)) == r


def test_combination_rank_unrank():
    for n in range(9):
        for k in range(n + 1):
            for r, c in enumerate(itertools.combinations(range(n), k)):
                assert combination_unrank(r, n, k) == list(c)
                assert combination_rank(c, n) == r


def test_validation():
    with pytest.raises(ValueError):
        permutation_unrank(10, 3)                        # rank >= 3!
    with pytest.raises(ValueError):
        combination_unrank(100, 5, 2)                    # rank >= C(5,2)
    with pytest.raises(ValueError):
        permutation_rank([0, 0, 1])                      # not a permutation
