"""Tests for Z-function, prefix function, periodicity, and Manacher, vs brute force."""

import random

from quantforge.string_periodicity import (
    z_function,
    prefix_function,
    smallest_period,
    is_periodic,
    borders,
    count_occurrences,
    manacher_longest_palindrome,
)


def _brute_z(s):
    n = len(s)
    z = [0] * n
    for i in range(1, n):
        k = 0
        while i + k < n and s[k] == s[i + k]:
            k += 1
        z[i] = k
    return z


def _brute_pi(s):
    n = len(s)
    pi = [0] * n
    for i in range(1, n):
        for L in range(i, 0, -1):
            if s[:L] == s[i - L + 1:i + 1]:
                pi[i] = L
                break
    return pi


def _brute_period(s):
    n = len(s)
    if n == 0:
        return 0
    for p in range(1, n + 1):
        if all(s[i] == s[i - p] for i in range(p, n)):
            return p
    return n


def _brute_occ(t, p):
    if p == "":
        return list(range(len(t) + 1))
    return [i for i in range(len(t) - len(p) + 1) if t[i:i + len(p)] == p]


def _brute_lps(s):
    best = ""
    for i in range(len(s)):
        for j in range(i + 1, len(s) + 1):
            sub = s[i:j]
            if sub == sub[::-1] and len(sub) > len(best):
                best = sub
    return best


def test_fuzz_z_and_prefix():
    rng = random.Random(121)
    for _ in range(4000):
        s = "".join(rng.choice("ab") for _ in range(rng.randint(0, 16)))
        assert z_function(s) == _brute_z(s)
        assert prefix_function(s) == _brute_pi(s)


def test_fuzz_period():
    rng = random.Random(122)
    for _ in range(4000):
        s = "".join(rng.choice("ab") for _ in range(rng.randint(0, 16)))
        assert smallest_period(s) == _brute_period(s)


def test_fuzz_palindrome_length_and_validity():
    rng = random.Random(123)
    for _ in range(4000):
        s = "".join(rng.choice("ab") for _ in range(rng.randint(0, 16)))
        lps = manacher_longest_palindrome(s)
        assert len(lps) == len(_brute_lps(s))
        assert lps == lps[::-1]
        assert lps in s


def test_fuzz_occurrences():
    rng = random.Random(124)
    for _ in range(4000):
        s = "".join(rng.choice("ab") for _ in range(rng.randint(0, 16)))
        p = "".join(rng.choice("ab") for _ in range(rng.randint(0, 4)))
        assert count_occurrences(s, p) == _brute_occ(s, p)


def test_prefix_function_known():
    assert prefix_function("abcabcd") == [0, 0, 0, 1, 2, 3, 0]
    assert prefix_function("aabaaab") == [0, 1, 0, 1, 2, 2, 3]


def test_smallest_period_and_is_periodic():
    assert smallest_period("abcabcabc") == 3
    assert is_periodic("abcabcabc") is True
    assert smallest_period("abcabca") == 3  # partial last block
    assert is_periodic("abcabca") is False
    assert smallest_period("abcd") == 4
    assert is_periodic("abcd") is False
    assert smallest_period("aaaa") == 1
    assert is_periodic("aaaa") is True


def test_borders():
    assert borders("abacaba") == [1, 3]  # "a", "aba"
    assert borders("aaaa") == [1, 2, 3]
    assert borders("abcd") == []


def test_manacher_known():
    assert manacher_longest_palindrome("babad") in ("bab", "aba")
    assert manacher_longest_palindrome("cbbd") == "bb"
    assert manacher_longest_palindrome("a") == "a"
    assert manacher_longest_palindrome("") == ""
    assert manacher_longest_palindrome("abccba") == "abccba"


def test_manacher_ties_earliest():
    # "abaxyzzyx" -> longest is "xyzzyx" (len 6); ensure a clean unique-length case
    assert manacher_longest_palindrome("forgeeksskeegfor") == "geeksskeeg"


def test_count_occurrences_overlap():
    assert count_occurrences("ababab", "ab") == [0, 2, 4]
    assert count_occurrences("aaa", "aa") == [0, 1]
    assert count_occurrences("abc", "xyz") == []


def test_empty_inputs():
    assert z_function("") == []
    assert prefix_function("") == []
    assert smallest_period("") == 0
    assert is_periodic("") is False
    assert borders("") == []
    assert manacher_longest_palindrome("") == ""


def test_empty_pattern_matches_everywhere():
    assert count_occurrences("abc", "") == [0, 1, 2, 3]


def test_pattern_longer_than_text():
    assert count_occurrences("ab", "abc") == []
