"""Tests for AhoCorasick multi-pattern matching.

Cross-checked against a brute-force per-pattern search that allows overlaps.
"""

import random

import pytest

from quantforge.aho_corasick import AhoCorasick


def _brute(patterns, text):
    out = []
    for p in set(patterns):
        if not p:
            continue
        start = 0
        while True:
            idx = text.find(p, start)
            if idx < 0:
                break
            out.append((idx + len(p) - 1, p))
            start = idx + 1  # allow overlaps
    return sorted(out)


def test_fuzz_matches_brute_force():
    rng = random.Random(11)
    alpha = "abc"
    for _ in range(2000):
        pats = [
            "".join(rng.choice(alpha) for _ in range(rng.randint(1, 4)))
            for _ in range(rng.randint(1, 5))
        ]
        text = "".join(rng.choice(alpha) for _ in range(rng.randint(0, 30)))
        ac = AhoCorasick(pats)
        assert sorted(ac.find_all(text)) == _brute(pats, text)


def test_classic_ushers_example():
    ac = AhoCorasick(["he", "she", "his", "hers"])
    res = sorted(ac.find_all("ushers"))
    assert (3, "she") in res
    assert (3, "he") in res
    assert (5, "hers") in res
    assert (0, "his") not in res


def test_overlapping_repeats():
    ac = AhoCorasick(["aa"])
    assert ac.find_all("aaaa") == [(1, "aa"), (2, "aa"), (3, "aa")]


def test_output_link_reports_shorter_pattern():
    # "he" is a suffix of "she"; both must be reported when "she" matches
    ac = AhoCorasick(["she", "he"])
    res = sorted(ac.find_all("she"))
    assert (2, "she") in res
    assert (2, "he") in res


def test_contains_any():
    ac = AhoCorasick(["cat", "dog"])
    assert ac.contains_any("the dog ran") is True
    assert ac.contains_any("nothing here") is False


def test_count_matches_counts_overlaps():
    ac = AhoCorasick(["cat", "dog"])
    assert ac.count_matches("catdogcat") == 3
    assert AhoCorasick(["aba"]).count_matches("ababa") == 2


def test_duplicate_and_empty_patterns_collapsed():
    ac = AhoCorasick(["ab", "ab", "", ""])
    assert ac.patterns == ["ab"]
    assert ac.find_all("abab") == [(1, "ab"), (3, "ab")]


def test_no_patterns():
    assert AhoCorasick([]).find_all("abc") == []
    assert AhoCorasick([]).contains_any("abc") is False
    assert AhoCorasick([]).count_matches("abc") == 0


def test_empty_text():
    assert AhoCorasick(["x"]).find_all("") == []
    assert AhoCorasick(["x"]).contains_any("") is False


def test_pattern_longer_than_text():
    assert AhoCorasick(["abcdef"]).find_all("abc") == []


def test_add_after_build_raises():
    ac = AhoCorasick(["a"])
    with pytest.raises(RuntimeError):
        ac.add("b")


def test_incremental_build():
    ac = AhoCorasick()
    ac.add("foo").add("bar")
    ac.build()
    res = sorted(ac.find_all("foobar"))
    assert res == [(2, "foo"), (5, "bar")]


def test_find_is_lazy_generator():
    ac = AhoCorasick(["ab"])
    gen = ac.find("abab")
    first = next(gen)
    assert first == (1, "ab")


def test_reuse_across_multiple_texts():
    ac = AhoCorasick(["ing", "ed"])
    assert ac.count_matches("running") == 1
    assert ac.count_matches("jumped") == 1
    assert ac.count_matches("singing") == 2


def test_patterns_property_returns_copy():
    ac = AhoCorasick(["a", "b"])
    p = ac.patterns
    p.append("c")
    assert "c" not in ac.patterns
