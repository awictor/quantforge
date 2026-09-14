"""Sequence alignment: Damerau-Levenshtein and Needleman-Wunsch."""

import random
from collections import deque

import pytest

from quantforge import damerau_levenshtein, needleman_wunsch
from quantforge import levenshtein


def test_transposition_is_one_edit():
    assert damerau_levenshtein("ca", "ac") == 1
    assert levenshtein("ca", "ac") == 2


def test_dl_at_most_levenshtein():
    rng = random.Random(1)
    for _ in range(2000):
        a = "".join(rng.choice("abc") for _ in range(rng.randint(0, 7)))
        b = "".join(rng.choice("abc") for _ in range(rng.randint(0, 7)))
        assert damerau_levenshtein(a, b) <= levenshtein(a, b)
        assert damerau_levenshtein(a, a) == 0


def test_dl_vs_brute_bfs():
    def brute(a, b, alpha="ab"):
        if a == b:
            return 0
        seen = {a}
        q = deque([(a, 0)])
        while q:
            s, d = q.popleft()
            if s == b:
                return d
            if d > 4:
                continue
            nb = set()
            for i in range(len(s) + 1):
                for c in alpha:
                    nb.add(s[:i] + c + s[i:])
            for i in range(len(s)):
                nb.add(s[:i] + s[i + 1:])
                for c in alpha:
                    nb.add(s[:i] + c + s[i + 1:])
            for i in range(len(s) - 1):
                nb.add(s[:i] + s[i + 1] + s[i] + s[i + 2:])
            for t in nb:
                if t not in seen:
                    seen.add(t)
                    q.append((t, d + 1))
        return None

    rng = random.Random(2)
    for _ in range(300):
        a = "".join(rng.choice("ab") for _ in range(rng.randint(0, 4)))
        b = "".join(rng.choice("ab") for _ in range(rng.randint(0, 4)))
        assert damerau_levenshtein(a, b) == brute(a, b)


def test_needleman_wunsch():
    def score(aa, bb, match=1, mismatch=-1, gap=-1):
        s = 0
        for x, y in zip(aa, bb):
            if x == "-" or y == "-":
                s += gap
            elif x == y:
                s += match
            else:
                s += mismatch
        return s

    rng = random.Random(3)
    for _ in range(1000):
        a = "".join(rng.choice("acgt") for _ in range(rng.randint(1, 10)))
        b = "".join(rng.choice("acgt") for _ in range(rng.randint(1, 10)))
        sc, aa, bb = needleman_wunsch(a, b)
        assert len(aa) == len(bb)
        assert score(aa, bb) == sc
        assert aa.replace("-", "") == a and bb.replace("-", "") == b
    assert needleman_wunsch("hello", "hello")[0] == 5
