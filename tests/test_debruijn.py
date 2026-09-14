"""Tests for de Bruijn sequences and Lyndon words, cross-checked against brute force."""

import itertools

import pytest

from quantforge.debruijn import de_bruijn_sequence, lyndon_words, is_lyndon


def _brute_lyndon(k, max_len):
    out = []
    for L in range(1, max_len + 1):
        for w in itertools.product(range(k), repeat=L):
            if is_lyndon(w):
                out.append(w)
    return sorted(out)


def test_de_bruijn_covers_every_tuple_once():
    for k in range(1, 5):
        for n in range(1, 6):
            if k ** n > 20000:
                continue
            seq = de_bruijn_sequence(k, n)
            assert len(seq) == k ** n
            L = len(seq)
            seen = {}
            for i in range(L):
                tup = tuple(seq[(i + j) % L] for j in range(n))
                seen[tup] = seen.get(tup, 0) + 1
            assert set(seen) == set(itertools.product(range(k), repeat=n))
            assert all(v == 1 for v in seen.values())


def test_lyndon_words_vs_brute():
    for k in range(1, 5):
        for max_len in range(0, 6):
            assert sorted(lyndon_words(k, max_len)) == _brute_lyndon(k, max_len)


def test_is_lyndon_explicit():
    assert is_lyndon([0, 0, 1]) is True
    assert is_lyndon([0, 1]) is True
    assert is_lyndon([0, 0, 1, 0, 1, 1]) is True
    assert is_lyndon([0, 1, 0, 1]) is False  # equals a rotation of itself (periodic)
    assert is_lyndon([1, 0]) is False
    assert is_lyndon([1, 1]) is False
    assert is_lyndon([]) is False
    assert is_lyndon([0]) is True  # single symbol is Lyndon


def test_de_bruijn_b2_3():
    s = de_bruijn_sequence(2, 3)
    assert len(s) == 8
    assert set(s) <= {0, 1}


def test_de_bruijn_b2_2():
    s = de_bruijn_sequence(2, 2)
    assert len(s) == 4
    L = 4
    tuples = {tuple(s[(i + j) % L] for j in range(2)) for i in range(L)}
    assert tuples == {(0, 0), (0, 1), (1, 0), (1, 1)}


def test_de_bruijn_unary():
    # k=1: only symbol 0, sequence is a single 0 repeated conceptually; length 1
    assert de_bruijn_sequence(1, 1) == [0]
    assert de_bruijn_sequence(1, 3) == [0]


def test_de_bruijn_large_alphabet():
    seq = de_bruijn_sequence(4, 2)
    assert len(seq) == 16
    L = 16
    tuples = {tuple(seq[(i + j) % L] for j in range(2)) for i in range(L)}
    assert tuples == set(itertools.product(range(4), repeat=2))


def test_lyndon_count_matches_necklace_formula():
    # number of Lyndon words of length exactly n over k symbols: (1/n) sum_{d|n} mu(d) k^(n/d)
    # for k=2, n=1..6: 2, 1, 2, 3, 6, 9
    counts = {}
    for w in lyndon_words(2, 6):
        counts[len(w)] = counts.get(len(w), 0) + 1
    assert [counts[n] for n in range(1, 7)] == [2, 1, 2, 3, 6, 9]


def test_invalid_arguments_raise():
    with pytest.raises(ValueError):
        de_bruijn_sequence(0, 2)
    with pytest.raises(ValueError):
        de_bruijn_sequence(2, 0)
    with pytest.raises(ValueError):
        lyndon_words(0, 3)
    with pytest.raises(ValueError):
        lyndon_words(2, -1)
