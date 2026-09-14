"""Tests for Trie (prefix tree), cross-checked against a plain set/list reference."""

import random

from quantforge.trie import Trie


def _ref_count(words, p):
    return sum(1 for w in words if w.startswith(p))


def _ref_keys(words, p):
    return sorted(w for w in words if w.startswith(p))


def _ref_longest(words, q):
    best = ""
    for w in words:
        if q.startswith(w) and len(w) > len(best):
            best = w
    return best


def test_fuzz_vs_reference_with_deletes():
    rng = random.Random(111)
    for _ in range(3000):
        words = set()
        for _ in range(rng.randint(0, 20)):
            words.add("".join(rng.choice("abc") for _ in range(rng.randint(0, 5))))
        t = Trie(words)
        assert len(t) == len(words)
        for _ in range(5):
            q = "".join(rng.choice("abc") for _ in range(rng.randint(0, 5)))
            assert t.contains(q) == (q in words)
            assert (q in t) == (q in words)
        for _ in range(5):
            p = "".join(rng.choice("abc") for _ in range(rng.randint(0, 4)))
            assert t.count_prefix(p) == _ref_count(words, p)
            assert t.starts_with(p) == (_ref_count(words, p) > 0)
            assert t.keys_with_prefix(p) == _ref_keys(words, p)
            q = "".join(rng.choice("abc") for _ in range(rng.randint(0, 6)))
            assert t.longest_prefix_of(q) == _ref_longest(words, q)
        assert t.keys() == sorted(words)
        # interleaved deletes
        wl = list(words)
        rng.shuffle(wl)
        for w in wl[: len(wl) // 2]:
            assert t.delete(w) is True
            words.discard(w)
        assert t.delete("ZZZ_absent") is False
        assert len(t) == len(words)
        assert t.keys() == sorted(words)


def test_explicit_prefix_queries():
    t = Trie(["cat", "car", "card", "dog"])
    assert t.keys_with_prefix("ca") == ["car", "card", "cat"]
    assert t.count_prefix("car") == 2
    assert t.starts_with("do") is True
    assert t.starts_with("z") is False
    assert t.contains("car") is True
    assert t.contains("ca") is False  # a prefix, not a stored key


def test_longest_prefix_of():
    t = Trie(["cat", "car", "card"])
    assert t.longest_prefix_of("cards") == "card"
    assert t.longest_prefix_of("cart") == "car"
    assert t.longest_prefix_of("ca") == ""  # "ca" is not stored
    assert t.longest_prefix_of("x") == ""


def test_longest_prefix_with_empty_key_stored():
    t = Trie(["", "ab"])
    assert t.longest_prefix_of("abc") == "ab"
    assert t.longest_prefix_of("zzz") == ""  # empty string is the best match


def test_empty_string_key():
    t = Trie([""])
    assert t.contains("") is True
    assert t.count_prefix("") == 1
    assert len(t) == 1
    assert t.keys() == [""]


def test_empty_trie():
    t = Trie()
    assert len(t) == 0
    assert t.keys() == []
    assert t.starts_with("a") is False
    assert t.count_prefix("") == 0
    assert t.longest_prefix_of("abc") == ""
    assert t.delete("x") is False


def test_double_insert():
    t = Trie()
    assert t.insert("a") is True
    assert t.insert("a") is False
    assert len(t) == 1


def test_delete_prunes_but_keeps_shared_prefix():
    t = Trie(["car", "card"])
    assert t.delete("card") is True
    assert t.contains("car") is True
    assert t.contains("card") is False
    assert t.count_prefix("car") == 1


def test_delete_absent_returns_false():
    t = Trie(["cat"])
    assert t.delete("dog") is False
    assert t.delete("ca") is False  # prefix but not a key
    assert len(t) == 1


def test_reinsert_after_delete():
    t = Trie(["cat"])
    t.delete("cat")
    assert not t.contains("cat")
    assert t.insert("cat") is True
    assert t.contains("cat")
    assert len(t) == 1


def test_count_prefix_empty_prefix_is_size():
    t = Trie(["a", "b", "ab", "abc"])
    assert t.count_prefix("") == 4
    assert t.count_prefix("ab") == 2


def test_keys_with_prefix_no_match():
    t = Trie(["cat", "dog"])
    assert t.keys_with_prefix("x") == []


def test_len_and_contains_dunder():
    t = Trie(["x", "y"])
    assert len(t) == 2
    assert "x" in t
    assert "z" not in t
