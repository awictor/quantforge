"""Fuzzy string similarity: Jaro, Jaro-Winkler, Dice, Jaccard."""

import random

import pytest

from quantforge import jaro, jaro_winkler, dice_coefficient, jaccard_similarity


def test_jaro_reference_values():
    assert abs(jaro("MARTHA", "MARHTA") - 0.944) < 1e-3
    assert abs(jaro("DWAYNE", "DUANE") - 0.822) < 1e-3
    assert abs(jaro("DIXON", "DICKSONX") - 0.767) < 1e-3


def test_jaro_winkler_reference_values():
    assert abs(jaro_winkler("MARTHA", "MARHTA") - 0.961) < 1e-3
    assert abs(jaro_winkler("DWAYNE", "DUANE") - 0.84) < 1e-3
    assert abs(jaro_winkler("DIXON", "DICKSONX") - 0.813) < 1e-3


def test_identity_and_empty():
    assert jaro("abc", "abc") == 1.0
    assert jaro_winkler("abc", "abc") == 1.0
    assert jaro("", "abc") == 0.0
    assert jaro("", "") == 1.0


def test_jaro_symmetric_and_bounded():
    rng = random.Random(1)
    for _ in range(500):
        a = "".join(rng.choice("abcde") for _ in range(rng.randint(0, 8)))
        b = "".join(rng.choice("abcde") for _ in range(rng.randint(0, 8)))
        assert abs(jaro(a, b) - jaro(b, a)) < 1e-12
        assert 0.0 <= jaro(a, b) <= 1.0 + 1e-9
        assert 0.0 <= jaro_winkler(a, b) <= 1.0 + 1e-9
        assert jaro_winkler(a, b) >= jaro(a, b) - 1e-12   # prefix boost never hurts


def test_dice_coefficient():
    assert abs(dice_coefficient("night", "nacht") - 0.25) < 1e-9
    assert dice_coefficient("abc", "abc") == 1.0
    assert dice_coefficient("ab", "xy") == 0.0


def test_jaccard_similarity():
    assert jaccard_similarity("the cat sat", "the dog sat") == 0.5
    assert abs(jaccard_similarity("abc", "bcd", tokenize=list) - 0.5) < 1e-9
    assert jaccard_similarity("", "") == 1.0
