"""Smith-Waterman local alignment."""

import random

from quantforge import smith_waterman


def _score(aa, bb, match=2, mismatch=-1, gap=-1):
    s = 0
    for x, y in zip(aa, bb):
        if x == "-" or y == "-":
            s += gap
        elif x == y:
            s += match
        else:
            s += mismatch
    return s


def test_embedded_substring():
    sc, aa, bb = smith_waterman("xxxHELLOyyy", "zzHELLOww")
    assert sc == 10
    assert aa.replace("-", "") == "HELLO" and bb.replace("-", "") == "HELLO"


def test_score_and_substring_property():
    rng = random.Random(1)
    for _ in range(1000):
        a = "".join(rng.choice("acgt") for _ in range(rng.randint(1, 12)))
        b = "".join(rng.choice("acgt") for _ in range(rng.randint(1, 12)))
        sc, aa, bb = smith_waterman(a, b)
        assert sc >= 0 and len(aa) == len(bb)
        if sc > 0:
            assert _score(aa, bb) == sc
        assert aa.replace("-", "") in a and bb.replace("-", "") in b


def test_no_common_and_identical():
    assert smith_waterman("aaa", "ttt") == (0, "", "")
    assert smith_waterman("hello", "hello")[0] == 10
