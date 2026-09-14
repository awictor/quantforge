"""String algorithms: edit distance, LCS, substring, KMP search."""

import random

import pytest

from quantforge import (
    levenshtein,
    hamming_distance,
    longest_common_subsequence,
    longest_common_substring,
    kmp_search,
)


def test_levenshtein_known():
    assert levenshtein("kitten", "sitting") == 3
    assert levenshtein("flaw", "lawn") == 2
    assert levenshtein("abc", "abc") == 0
    assert levenshtein("", "abc") == 3


def test_levenshtein_matches_brute_and_symmetric():
    def brute(a, b):
        if not a:
            return len(b)
        if not b:
            return len(a)
        if a[0] == b[0]:
            return brute(a[1:], b[1:])
        return 1 + min(brute(a[1:], b), brute(a, b[1:]), brute(a[1:], b[1:]))

    rng = random.Random(1)
    for _ in range(300):
        a = "".join(rng.choice("abc") for _ in range(rng.randint(0, 6)))
        b = "".join(rng.choice("abc") for _ in range(rng.randint(0, 6)))
        assert levenshtein(a, b) == brute(a, b) == levenshtein(b, a)


def test_lcs():
    assert len(longest_common_subsequence("ABCBDAB", "BDCAB")) == 4

    def is_subseq(s, t):
        it = iter(t)
        return all(c in it for c in s)

    rng = random.Random(2)
    for _ in range(200):
        a = "".join(rng.choice("abc") for _ in range(rng.randint(0, 8)))
        b = "".join(rng.choice("abc") for _ in range(rng.randint(0, 8)))
        l = longest_common_subsequence(a, b)
        assert is_subseq(l, a) and is_subseq(l, b)


def test_longest_common_substring():
    assert longest_common_substring("abcdxyz", "xyzabcd") == "abcd"
    assert longest_common_substring("abc", "xyz") == ""


def test_kmp_matches_naive():
    def naive(text, pat):
        return [i for i in range(len(text) - len(pat) + 1) if text[i:i + len(pat)] == pat]

    rng = random.Random(3)
    for _ in range(300):
        text = "".join(rng.choice("ab") for _ in range(rng.randint(0, 20)))
        pat = "".join(rng.choice("ab") for _ in range(rng.randint(1, 4)))
        assert kmp_search(text, pat) == naive(text, pat)
    assert kmp_search("aaaa", "aa") == [0, 1, 2]        # overlapping
    assert kmp_search("ab", "") == [0, 1, 2]            # empty pattern


def test_hamming():
    assert hamming_distance("karolin", "kathrin") == 3
    assert hamming_distance("abc", "abc") == 0
    with pytest.raises(ValueError):
        hamming_distance("ab", "abc")
