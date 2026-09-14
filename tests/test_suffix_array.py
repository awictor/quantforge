"""Tests for suffix array / LCP array and their substring queries.

Cross-checked against brute-force references (sorting suffixes, naive LCP, set of
substrings, overlapping-occurrence counting).
"""

import random

from quantforge.suffix_array import (
    suffix_array,
    rank_array,
    lcp_array,
    substring_search,
    longest_repeated_substring,
    count_distinct_substrings,
)


def _brute_sa(t):
    return sorted(range(len(t)), key=lambda i: t[i:])


def _brute_lcp(t, sa):
    n = len(t)
    out = [0] * n
    for r in range(1, n):
        a, b = t[sa[r - 1]:], t[sa[r]:]
        h = 0
        while h < len(a) and h < len(b) and a[h] == b[h]:
            h += 1
        out[r] = h
    return out


def _brute_search(t, p):
    if p == "":
        return list(range(len(t) + 1)) if t else [0]
    return [i for i in range(len(t) - len(p) + 1) if t[i:i + len(p)] == p]


def _occ_overlap(t, s):
    c = 0
    i = t.find(s)
    while i >= 0:
        c += 1
        i = t.find(s, i + 1)
    return c


def _brute_lrs_len(t):
    best = 0
    for i in range(len(t)):
        for j in range(i + 1, len(t) + 1):
            s = t[i:j]
            if len(s) > best and _occ_overlap(t, s) >= 2:
                best = len(s)
    return best


def _brute_distinct(t):
    return len({t[i:j] for i in range(len(t)) for j in range(i + 1, len(t) + 1)})


def test_suffix_array_fuzz_vs_brute():
    rng = random.Random(21)
    for _ in range(3000):
        t = "".join(rng.choice("ab") for _ in range(rng.randint(0, 14)))
        assert suffix_array(t) == _brute_sa(t)


def test_lcp_array_fuzz_vs_brute():
    rng = random.Random(22)
    for _ in range(2000):
        t = "".join(rng.choice("abc") for _ in range(rng.randint(0, 16)))
        sa = suffix_array(t)
        assert lcp_array(t, sa) == _brute_lcp(t, sa)


def test_substring_search_fuzz_vs_brute():
    rng = random.Random(23)
    for _ in range(2000):
        t = "".join(rng.choice("ab") for _ in range(rng.randint(0, 14)))
        sa = suffix_array(t)
        for _ in range(3):
            p = "".join(rng.choice("ab") for _ in range(rng.randint(0, 4)))
            assert substring_search(t, p, sa) == _brute_search(t, p)


def test_distinct_substrings_fuzz_vs_brute():
    rng = random.Random(24)
    for _ in range(2000):
        t = "".join(rng.choice("abc") for _ in range(rng.randint(0, 14)))
        assert count_distinct_substrings(t) == _brute_distinct(t)


def test_longest_repeated_length_matches_brute():
    rng = random.Random(25)
    for _ in range(2000):
        t = "".join(rng.choice("ab") for _ in range(rng.randint(0, 16)))
        lrs = longest_repeated_substring(t)
        assert len(lrs) == _brute_lrs_len(t)
        if lrs:
            assert _occ_overlap(t, lrs) >= 2


def test_banana_suffix_array():
    t = "banana"
    assert suffix_array(t) == [5, 3, 1, 0, 4, 2]


def test_banana_lcp():
    t = "banana"
    sa = suffix_array(t)
    # suffixes sorted: a, ana, anana, banana, na, nana
    assert lcp_array(t, sa) == [0, 1, 3, 0, 0, 2]


def test_banana_longest_repeated():
    assert longest_repeated_substring("banana") == "ana"


def test_banana_distinct_count():
    assert count_distinct_substrings("banana") == 15


def test_banana_search():
    assert substring_search("banana", "ana") == [1, 3]
    assert substring_search("banana", "na") == [2, 4]
    assert substring_search("banana", "xyz") == []


def test_rank_array_is_inverse():
    t = "mississippi"
    sa = suffix_array(t)
    rank = rank_array(sa)
    for r, i in enumerate(sa):
        assert rank[i] == r


def test_empty_and_singleton():
    assert suffix_array("") == []
    assert lcp_array("") == []
    assert suffix_array("x") == [0]
    assert lcp_array("x") == [0]
    assert longest_repeated_substring("") == ""
    assert longest_repeated_substring("x") == ""
    assert count_distinct_substrings("") == 0


def test_no_repeat_returns_empty():
    assert longest_repeated_substring("abcd") == ""


def test_all_same_character():
    t = "aaaa"
    assert suffix_array(t) == [3, 2, 1, 0]
    assert longest_repeated_substring(t) == "aaa"
    assert count_distinct_substrings(t) == 4  # a, aa, aaa, aaaa


def test_empty_pattern_search():
    assert substring_search("abc", "") == [0, 1, 2, 3]
    assert substring_search("", "") == [0]


def test_search_on_empty_text():
    assert substring_search("", "a") == []
